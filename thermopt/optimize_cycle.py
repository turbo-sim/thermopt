# Standard library
import os
import time
import copy
import pickle
import datetime
import threading

# Third-party
import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import jaxprop as cpx
import pysolver_view as psv

# Local modules
from . import cycles
from . import utilities as utils
from .config import read_configuration_file
from .cycles.state_requirements import validate_cycle_state_requirements
from .cycles.setup import (
    prepare_cycle_models, prepare_cycle_references, calculate_reference_points,
)

COLORS_MATLAB = utils.COLORS_MATLAB
LABEL_MAPPING = {
    "entropy": "Entropy [J/kg/K]",
    "temperature": "Temperature [K]",
    "pressure": "Pressure [Pa]",
    "enthalpy": "Enthalpy [J/kg]",
    "density": "Density [kg/m$^3$]",
    "speed_of_sound": "Speed of sound [m/s]",
    "compressibility_factor": "Compressibility factor [-]",
    "heat": "Heat flow rate [W]",
    "heat_flow": "Heat flow rate [W]",
}

# Cycle configurations available
CYCLE_TOPOLOGIES = {
    "simple": cycles.cycle_power_simple.evaluate_cycle,
    "power_simple": cycles.cycle_power_simple.evaluate_cycle,
    "recuperated": cycles.cycle_power_recuperated.evaluate_cycle,
    "power_recuperated": cycles.cycle_power_recuperated.evaluate_cycle,
    "power_split_compression": cycles.cycle_power_split_compression.evaluate_cycle,
    "split_compression": cycles.cycle_power_split_compression.evaluate_cycle,
    "recompression": cycles.cycle_power_split_compression.evaluate_cycle,
    "refrigeration_simple": cycles.cycle_refrigeration_simple.evaluate_cycle,
    "refrigeration_recuperated": cycles.cycle_refrigeration_recuperated.evaluate_cycle,
    "PTES_recuperated": cycles.cycle_PTES_recuperated.evaluate_cycle,
    "PTES_recuperated_turbo": cycles.cycle_PTES_recuperated_turbo.evaluate_cycle,
}

GRAPHICS_PLACEHOLDER = {
    "process_lines": {},
    "state_points": {},
    "pinch_point_lines": {},
}


class ThermodynamicCycleOptimization:
    def __init__(self, config_file, out_dir=None):
        """
        Initializes the optimization manager with a configuration file.
        Parameters:
            config_file (str): The path to the YAML configuration file.
        """

        # Create output directory
        if out_dir is None:
            timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
            self.out_dir = f"results/case_{timestamp}"
        else:
            self.out_dir = out_dir

        self.optimization_dir = os.path.join(self.out_dir, "optimization")
        os.makedirs(self.optimization_dir, exist_ok=True)

        # Read configuration file
        self.config = self.read_config(config_file)

    def read_config(self, config_file):
        """
        Loads configuration from a YAML file.
        """
        config = read_configuration_file(config_file)
        self.load_config(config)
        return self.config

    def load_config(self, config_dict):
        """
        Load a new configuration and update the problem and solver objects.

        Reuse this manager's problem and unchanged native models; reconstruct
        the solver so current constraint dimensions apply. Configuration is
        copied. If preparation, initial evaluation or solver setup fails, keep
        the previous configured problem, solver and numerical results.

        Parameters:
            config_obj (dict): A dictionary-like configuration object.
        """
        # Configuration and solver changes are orchestration transactions. A
        # failed reload must leave the last evaluated problem usable. Native
        # flash cursors are transient; immutable returned states remain valid.
        previous_config = getattr(self, "config", None)
        previous_problem = getattr(self, "problem", None)
        previous_solver = getattr(self, "solver", None)
        previous_state = None if previous_problem is None else previous_problem.__dict__.copy()
        try:
            self.config = copy.deepcopy(config_dict)
            self.problem = self.setup_problem()
            self.solver = self.setup_solver()
        except Exception:
            if previous_state is not None:
                previous_problem.__dict__.clear()
                previous_problem.__dict__.update(previous_state)
            self.config, self.problem, self.solver = previous_config, previous_problem, previous_solver
            raise
        return self.config

    def setup_problem(self):
        """
        Sets up the ThermodynamicCycleProblem based on the loaded configuration.
        """
        # Keep this manager's configured problem and its native models. Solver
        # objects are rebuilt so changed constraint dimensions are respected.
        if getattr(self, "problem", None) is None:
            self.problem = ThermodynamicCycleProblem(
                self.config["problem_formulation"], out_dir=self.out_dir
            )
        else:
            self.problem.update_problem(self.config["problem_formulation"])
        self.problem.fitness(self.problem.x0)
        return self.problem

    def setup_solver(self):
        """
        Configures and returns the optimization solver.
        """
        solver_options = {
            k: v for k, v in self.config["solver_options"].items() if k != "callbacks"
        }
        self.solver = psv.OptimizationSolver(
            self.problem,
            **solver_options,  # Pass all options except "callbacks"
            callback_functions=None,
            plot_scale_constraints="log",
            tolerance_check_cache=1e-10
        )
        return self.solver

    def set_config_value(self, path, value, reload=True):
        """
        Updates a value in the nested config dictionary using a dot-separated path.
        Optionally refreshes problem and solver.

        Example:
            set_config_value("problem_formulation.fixed_parameters.expander.efficiency", 0.8)
        """
        keys = path.split(".")
        candidate = copy.deepcopy(self.config)
        cfg = candidate
        for key in keys[:-1]:
            cfg = cfg[key]
        cfg[keys[-1]] = value

        if reload:
            self.load_config(candidate)
        else:
            self.config = candidate

    def set_config_values(self, updates: dict):
        """
        Updates multiple config values in one call, triggering a single reload.
        """
        candidate = copy.deepcopy(self.config)
        for path, value in updates.items():
            keys = path.split(".")
            target = candidate
            for key in keys[:-1]:
                target = target[key]
            target[keys[-1]] = value
        self.load_config(candidate)

    def set_constraint(self, variable, type=None, value=None, normalize=None):
        """
        Adds or updates a constraint for the given variable name.

        Parameters:
            variable (str): Full variable path, e.g., "$components.cooler_charge.temperature_difference"
            type (str, optional): Constraint type, e.g., ">", "<", "="
            value (float, optional): Target value
            normalize (bool or float, optional): Normalization value or flag
        """
        candidate = copy.deepcopy(self.config)
        constraints = candidate["problem_formulation"].get("constraints", [])
        if not isinstance(constraints, list):
            constraints = list(constraints)
        
        # Find existing constraint
        for c in constraints:
            if c.get("variable") == variable:
                if type is not None:
                    c["type"] = type
                if value is not None:
                    c["value"] = value
                if normalize is not None:
                    c["normalize"] = normalize
                break
        else:
            # Add new constraint
            new_constraint = {"variable": variable}
            if type is not None:
                new_constraint["type"] = type
            if value is not None:
                new_constraint["value"] = value
            if normalize is not None:
                new_constraint["normalize"] = normalize
            constraints.append(new_constraint)

        candidate["problem_formulation"]["constraints"] = constraints
        self.load_config(candidate)


    def run_optimization(self, x0=None):
        """
        Executes the optimization process.
        """
        callback_flags = self.config["solver_options"].get("callbacks", {})
        self.solver.callback_functions = []

        callback_registry = {
            "plot_cycle": self.plot_cycle_callback,
            "save_plot": self.save_plot_callback,
            "save_config": self.save_config_callback,
            "save_report": self.save_report_callback,
            "plot_convergence": self._init_convergence_callback,
        }

        # Plot convergence callback is treated as a special case
        for key, func in callback_registry.items():
            if callback_flags.get(key, False):
                func() if key == "plot_convergence" else self.solver.callback_functions.append(func)

        # Run the optimization for the specified or default initial guess
        if x0 is None:
            self.solver.solve(self.problem.x0)
        else:
            self.solver.solve(x0)

    def save_results(self):
        """
        Saves the results of the optimization, including configurations and output files.
        """
        filename = os.path.join(self.out_dir, "optimal_solution")
        self.problem.save_data_to_excel(filename=filename + ".xlsx")
        print("hey")
        self.problem.save_current_configuration(filename=filename + ".yaml")
        self.print_convergence_history(savefile=True)
        self.print_optimization_report(savefile=True)
        self.plot_convergence_history(savefile=True, showfig=False)

        # Plot final solution
        self.problem.plot_cycle()
        self.problem.figure.tight_layout(pad=1)
        self.problem.figure.suptitle(None)
        filename = os.path.join(self.out_dir, "optimal_solution.png")
        self.problem.figure.savefig(filename, dpi=500)
        plt.close(self.problem.figure)
        # self.save_solver_pickle()


    def save_solver_pickle(self):
        """
        Sanitize solver (including deeply nested problem) and save it to a pickle file.
        """
        filename = "optimization_solver"
        # utils.dump_object_structure(self.solver, log_file="testing.txt")
        utils.save_to_pickle(self.solver, filename=filename, out_dir=self.out_dir, timestamp=False)


    def plot_convergence_history(self, savefile=False, showfig=True):
        filename = "convergence_history"
        self.solver.plot_convergence_history(
            savefile=savefile,
            filename=filename,
            output_dir=self.out_dir,
            showfig=showfig,
        )

    def print_convergence_history(self, savefile=False):
        filename = "convergence_history.txt"
        self.solver.print_convergence_history(
            savefile=savefile, filename=filename, output_dir=self.out_dir, to_console=False
        )

    def print_optimization_report(self, savefile=False):
        filename = "optimization_report.txt"
        self.solver.print_optimization_report(
            savefile=savefile, filename=filename, output_dir=self.out_dir,
            include_design_variables=True,
            include_constraints=True,
            include_kkt_conditions=True,
            to_console=False,
        )

    def create_animation(self, format="both", fps=1):
        """
        Creates an animation from optimization history.

        Parameters
        ----------
        format : str, optional
            Format of animation ("gif", "mp4", or "both"), default is "both".
        duration : float, optional
            Duration of each frame in GIF (default is 0.5 sec).
        fps : int, optional
            Frames per second for MP4 (default is 10).
        """
        image_folder = os.path.join(self.out_dir, "optimization")
        gif_file = os.path.join(self.out_dir, "optimization_animation.gif")
        mp4_file = os.path.join(self.out_dir, "optimization_animation.mp4")

        # Check if the folder contains images
        image_files = sorted(f for f in os.listdir(image_folder) if f.endswith(".png"))
        if not image_files:
            print("No images found for animation. Skipping animation creation.")
            return

        # Call the utility functions
        if format in ["gif", "both"]:
            utils.create_gif(image_folder, gif_file, duration=len(image_files) / fps)
            print(f"GIF saved at {gif_file}")

        if format in ["mp4", "both"]:
            utils.create_mp4(image_folder, mp4_file, fps=1)
            print(f"MP4 saved at {mp4_file}")

    # ------------------------------------------------------------------------- #
    # -------------------------- Callback functions --------------------------- #
    # ------------------------------------------------------------------------- #
    def save_config_callback(self, x, iter):
        """
        A callback function to save the current configuration during optimization iterations.

        Parameters:
        - x : The current solution vector from the optimizer.
        - iter : The current optimization iteration count.

        This function acts as a bridge between the optimizer callback requirements and the
        existing `save_current_configuration` function.
        """
        # Call the existing function to save the configuration
        filename = os.path.join(self.optimization_dir, f"iteration_{iter:03d}.yaml")
        self.problem.save_current_configuration(filename)

    def save_report_callback(self, x, iter):
        """
        A callback function to save the current configuration during optimization iterations.

        Parameters:
        - x : The current solution vector from the optimizer.
        - iter : The current optimization iteration count.

        This function acts as a bridge between the optimizer callback requirements and the
        existing `save_current_configuration` function.
        """
        self.solver.print_optimization_report(
            self.problem.scale_normalized_to_physical(x),
            include_kkt_conditions=True,
            savefile=True,
            filename=f"iteration_{iter:03d}.txt",
            output_dir=self.optimization_dir,
            to_console=False
        )

    def plot_cycle_callback(self, x, iter):
        """
        Plot the thermodynamic cycle during optimization.
        """
        self.problem.plot_cycle()
        self.problem.figure.suptitle(
            f"Optimization iteration: {iter:03d}", fontsize=14, y=0.95
        )
        self.problem.figure.tight_layout(pad=1)

    def save_plot_callback(self, x, iter):
        """
        Save the thermodynamic cycle figure during optimization.
        """
        filename = os.path.join(self.optimization_dir, f"iteration_{iter:03d}.png")
        self.problem.figure.savefig(filename, dpi=500)

    def _init_convergence_callback(self):
        self.solver._plot_convergence_callback([], [], initialize=True)
        self.solver.callback_functions.append(self.solver._plot_convergence_callback)



class ThermodynamicCycleProblem(psv.OptimizationProblem):
    """
    A class to represent a thermodynamic cycle optimization problem.

    This class provides functionalities to load and update the configuration for
    the thermodynamic cycle, and to perform interactive plotting based on the current
    configuration.

    Attributes
    ----------
    plot_initialized : bool
        Flag to indicate if the plot has been initialized.
    constraints : dict
        Dictionary holding the constraints of the problem.
    fixed_parameters : dict
        Dictionary holding the fixed parameters of the problem.
    design_variables : dict
        Dictionary holding the current values of the design variables.
    lower_bounds : dict
        Dictionary holding the lower bounds of the design variables.
    upper_bounds : dict
        Dictionary holding the upper bounds of the design variables.
    keys : list
        List of keys (names) of the design variables.
    x0 : np.ndarray
        Initial guess for the optimization problem.

    Methods
    -------
    update_config(configuration):
        Update the problem's configuration based on the provided dictionary.
    load_config_from_file(configuration_file):
        Load and update the problem's configuration from a specified file.
    plot_cycle_interactive(configuration_file, update_interval=0.20):
        Perform interactive plotting, updating the plot based on the configuration file.
    """

    def __init__(self, configuration, out_dir=None):
        """
        Constructs all the necessary attributes for the ThermodynamicCycleProblem object.

        Parameters
        ----------
        config : dict
            Dictionary containing the configuration for the thermodynamic cycle problem.
        """

        # TODO improve output directory functionality. Is it needed?

        # As the first step, process the fixed parameters for dynamic calculations
        self.fixed_parameters = configuration["fixed_parameters"]

        # Initialize variables
        self.figure = None
        self.figure_TQ = None
        self.configuration = copy.deepcopy(configuration)
        self.graphics = copy.deepcopy(GRAPHICS_PLACEHOLDER)

        # Update problem based on current configuration
        self.update_problem(self.configuration)

        # Define filename with unique date-time identifier
        self.out_dir = out_dir
        if self.out_dir is None:
            current_time = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
            self.out_dir = f"results/case_{current_time}"

        # Create a directory to save simulation results
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)

    def load_configuration_file(self, config_file):
        """
        Load and update the problem's configuration from a specified file.

        Useful to plot the cycle according to the latest version of the configuration
        file in real time (interactive initial guess generation)

        Parameters
        ----------
        config_file : str
            Path to the configuration file.
        """
        config = read_configuration_file(config_file)
        self.update_problem(config["problem_formulation"])

    def update_problem(self, configuration):
        """
        Update the problem's configuration based on the provided dictionary.

        Prepare a private candidate and apply it only after setup succeeds.
        Unchanged stream models and reference points are reused within this
        problem; source-provenance checks occur here, never during fitness.
        A successful update invalidates previously evaluated cycle results.

        Parameters
        ----------
        config : dict
            Dictionary containing the new configuration for the thermodynamic cycle problem.
        """
        # Build a candidate using explicit setup functions; publish it only
        # after model validation, reference preparation and rendering succeed.
        conf = copy.deepcopy(configuration)
        parameters = conf["fixed_parameters"]
        fluids, settings, sources = prepare_cycle_models(
            parameters, getattr(self, "fluids", None),
            getattr(self, "_fluid_settings", None), getattr(self, "_fluid_sources", None),
        )
        validate_cycle_state_requirements(
            conf["cycle_topology"], parameters, fluids,
            conf["constraints"], conf["objective_function"],
        )
        # Reference-only setup queries can lie outside a table's declared domain.
        # Evaluate those through the same source EOS, never through interpolation.
        reference, reference_key, points, params = prepare_cycle_references(
            parameters, fluids["working_fluid"], sources["working_fluid"],
            getattr(self, "reference_fluid", None), getattr(self, "_reference_key", None),
            getattr(self, "_reference_points", None),
        )
        names = list(conf["design_variables"])
        lower = {k: v["min"] for k, v in conf["design_variables"].items()}
        upper = {k: v["max"] for k, v in conf["design_variables"].items()}
        initial = {k: v["value"] for k, v in conf["design_variables"].items()}

        # Evaluate symbolic expressions using data in "params" dict
        lb = [utils.render_and_evaluate(lower[k], params) for k in names]
        ub = [utils.render_and_evaluate(upper[k], params) for k in names]
        x0 = np.array([utils.render_and_evaluate(initial[k], params) for k in names])
        plots = conf["plot_settings"]
        plot_signature = (
            conf["cycle_topology"], bool(plots.get("pinch_point_diagram", False)),
            tuple(tuple(diagram.get(key) for key in ("x_prop", "y_prop", "x_scale", "y_scale"))
                  for diagram in plots["diagrams"]), sources["working_fluid"],
        )
        self.__dict__.update(
            configuration=conf, cycle_topology=conf["cycle_topology"], plot_settings=plots,
            constraints=conf["constraints"], fixed_parameters=parameters,
            objective_function=conf["objective_function"], variable_names=names,
            lb_dict=lower, ub_dict=upper, x0_dict=initial, lb=lb, ub=ub, x0=x0,
            fluids=fluids, fluid=fluids["working_fluid"], reference_fluid=reference, params=params,
            _fluid_settings=settings, _fluid_sources=sources, _reference_key=reference_key,
            _reference_points=points, _plot_signature=plot_signature,
        )
        # Results describe a particular evaluated configuration. A successful
        # reload requires a new fitness evaluation before those results exist.
        for name in ("cycle_data", "f", "c_eq", "c_ineq", "constraint_data"):
            self.__dict__.pop(name, None)

    def _calculate_special_points(self):
        """
        Calculates and stores key thermodynamic states for the fluid used in the cycle.

        This method computes two specific thermodynamic states:
        1. Saturated liquid state at the heat sink inlet temperature.
        2. Dilute gas state at the heat source inlet temperature.

        If the heat sink inlet temperature is below the fluid's critical temperature,
        the saturated liquid state is calculated at this temperature. Above the
        critical temperature, the established initialization uses the critical-
        density isochore, rather than a heat-capacity-maximum pseudocritical curve.

        The dilute gas state is calculated at the heat source inlet temperature and at
        very low pressure (slightly above the liquid lower-endpoint pressure,
        which is the triple pressure for a pure fluid).

        The calculated states, along with the fluid's critical and triple point properties,
        are stored separately in `params['working_fluid']` for bound rendering.
        Fluid construction recipes in `fixed_parameters` are not overwritten.
        Pseudo-pure endpoint records retain their source interpretation.

        These states are intended to define the bounds of the design variables specified
        in the YAML configuration file, for example:

            compressor_inlet_enthalpy:
            min: 0.9*$working_fluid.liquid_at_ambient_temperature.h
            max: 2.0*$working_fluid.liquid_at_ambient_temperature.h
            value: 1.1*$working_fluid.liquid_at_ambient_temperature.h
            expander_inlet_pressure:
            min: 0.75*$working_fluid.critical_point.p
            max: 5.00*$working_fluid.critical_point.p
            value: 2.0*$working_fluid.critical_point.p
            expander_inlet_enthalpy:
            min: 1.10*$working_fluid.critical_point.h
            max: 1.00*$working_fluid.gas_at_maximum_temperature.h
            value: 0.90*$working_fluid.gas_at_maximum_temperature.h

        """
        # Explicit recalculation remains available for callers of this helper;
        # ordinary reloads reuse source points through prepare_cycle_references.
        self._reference_points = calculate_reference_points(self.fixed_parameters, self.reference_fluid)
        self.params = copy.deepcopy(self.fixed_parameters)
        self.params["working_fluid"] = dict(self._reference_points)

    def fitness(self, x):
        """
        Evaluate optimization problem
        """

        # Link variable names and values
        self.x0_dict = dict(zip(self.variable_names, x))

        # Update configuration with the current values of x
        for k, v in self.x0_dict.items():
            if k in self.configuration["design_variables"]:
                self.configuration["design_variables"][k]["value"] = v
            else:
                # Optionally handle the error or log a warning if the key does not exist
                raise KeyError(f"{k} is not a recognized design variable.")

        # Evaluate thermodynamic cycle
        if self.cycle_topology in CYCLE_TOPOLOGIES.keys():
            self.cycle_data = CYCLE_TOPOLOGIES[self.cycle_topology](
                self.x0_dict,
                self.fixed_parameters,
                self.constraints,
                self.objective_function,
                fluids=self.fluids,
            )
        else:
            options = ", ".join(f"'{k}'" for k in CYCLE_TOPOLOGIES.keys())
            raise ValueError(
                f"Invalid cycle topology: '{self.cycle_topology}'. Available options: {options}"
            )

        # Define objective function and constraints
        self.f = self.cycle_data["objective_function"]
        self.c_eq = self.cycle_data["equality_constraints"]
        self.c_ineq = self.cycle_data["inequality_constraints"]
        self.constraint_data = self.cycle_data["constraints_report"]
        # self.constraint_data_ = self.cycle_data["constraints_report"]

        # Combine objective function and constraints
        out = psv.combine_objective_and_constraints(self.f, self.c_eq, self.c_ineq)

        return out

    def get_bounds(self):
        return self.lb, self.ub

    def get_nec(self):
        return psv.count_constraints(self.c_eq)

    def get_nic(self):
        return psv.count_constraints(self.c_ineq)

    def save_current_configuration(self, filename):
        """Save the current configuration to a YAML file."""
        config_data = {k: v for k, v in self.configuration.items()}
        config_data = utils.convert_numpy_to_python(config_data, precision=12)
        with open(filename, "w") as file:
            yaml.dump(config_data, file, default_flow_style=False, sort_keys=False)

    def plot_cycle(self):
        """
        Plots or updates the thermodynamic cycle diagrams based on current settings,
        including the option to include a pinch point diagram.

        This function is capable of both creating new cycle diagrams and updating existing ones.
        It's particularly useful in dynamic scenarios such as during optimization steps,
        where the plot needs to be refreshed continually with new data.
        The method also supports real-time updates based on the latest configuration settings.

        The function first determines the number of subplots required based on the cycle diagrams
        specified in the 'plot_settings' attribute and whether a pinch point diagram is included.
        It then either initializes a new figure and axes or updates existing ones.
        Each thermodynamic diagram (phase diagram and cycle components) and the optional pinch point
        diagram are plotted or updated accordingly.

        The method ensures that the plot reflects the current state of the cycle, including any
        changes during optimization or adjustments to configuration settings.

        The data of the plots can be updated interactively, but the subplot objects created can
        only be specified upon class initialization. For example, it is not possible to switch
        on and off the pinch point diagram or to add new thermodynamic diagrams. The class would
        have to be re-initialized upon those scenarios. The reason for this choice is that re-creating
        a figure would be too time consuming and not practical for the real-time updating of the plots

        """

        if "cycle_data" not in self.__dict__:
            raise RuntimeError("Evaluate fitness after configuring or updating the problem before plotting cycle results.")

        # Rebuild source backgrounds/layout only when the next plot is drawn.
        # A failed manager reload therefore cannot destroy the last valid plot.
        if getattr(self, "_rendered_plot_signature", None) != self._plot_signature:
            for figure in (self.figure, self.figure_TQ):
                if figure is not None:
                    plt.close(figure)
            self.figure = self.figure_TQ = None
            self.graphics = copy.deepcopy(GRAPHICS_PLACEHOLDER)

        # Determine the number of subplots
        include_pinch_diagram = self.plot_settings.get("pinch_point_diagram", False)
        ncols = len(self.plot_settings["diagrams"]) + int(include_pinch_diagram)
        nrows = 1

        # Initialize the figure and axes
        if not (self.figure and plt.fignum_exists(self.figure.number)):

            # Reset the graphics objects if the figure was closed
            self.graphics = copy.deepcopy(GRAPHICS_PLACEHOLDER)
            self.figure, self.axes = plt.subplots(
                nrows, ncols, figsize=(5.2 * ncols, 4.8)
            )
            self.axes = utils.ensure_iterable(self.axes)

        # Plot or update each thermodynamic_diagram
        for i, ax in enumerate(self.axes[:-1] if include_pinch_diagram else self.axes):
            plot_config = self.plot_settings["diagrams"][i]
            plot_config = self._get_diagram_default_settings(plot_config)
            self._plot_thermodynamic_diagram(ax, self.cycle_data, plot_config, i)

        # Plot or update the pinch point diagram
        if include_pinch_diagram:
            self._plot_pinch_point_diagram(self.axes[-1], ncols - 1)

        # Adjust layout and refresh plot
        self.figure.tight_layout(pad=1)
        plt.draw()
        plt.pause(0.01)
        self._rendered_plot_signature = self._plot_signature

    def plot_cycle_realtime(
        self, configuration_file, update_interval=0.1, write_report=False
    ):
        """
        Perform interactive plotting, updating the plot based on the configuration file.

        Parameters
        ----------
        config_file : str
            Path to the configuration file.
        update_interval : float, optional
            Time interval in seconds between plot updates (default is 0.1 seconds).

        """

        def wait_for_input():
            input()
            self.enter_pressed = True

        # Initialize secondary thread to wait for user input
        self.enter_pressed = False
        self.input_thread = threading.Thread(target=wait_for_input)
        self.input_thread.daemon = True
        self.input_thread.start()

        # Print instructions message
        print("-" * 80)
        print(" Creating thermodynamic cycle interactive plot")
        print("-" * 80)
        print(f" The current configuration file is: '{configuration_file}'")
        print(" Modify the configuration file and save it to update the plot")
        print(" Try to find a good initial guess for the optimization.")
        print(" Using a feasible initial guess improves optimization convergence.")
        print(" Press 'enter' to continue.")
        print("-" * 80)

        # Update the plot until termination signal
        while not self.enter_pressed:
            # Read the configuration file
            self.load_configuration_file(configuration_file)
            self.fitness(self.x0)
            self.plot_cycle()
            self.figure.suptitle(
                f"Iterative thermodynamic cycle configuration", fontsize=14, y=0.95
            )

            # Write optimization report to file
            # if write_report:
            #     report = self.make_optimization_report(self.x0)
            #     filename="initial_guess_report.txt"
            #     fullfile = os.path.join(self.out_dir, filename)
            #     with open(fullfile, "w") as f:
            #         f.write(report)

            # Wait for the specified interval before updating again
            time.sleep(update_interval)

            # Exit interactive plotting when the user closes the figure
            if not plt.fignum_exists(self.figure.number):
                break

            # Exit interactive plotting when the user presses enter
            if self.enter_pressed:
                plt.close(self.figure)
                break

    def _plot_thermodynamic_diagram(self, ax, cycle_data, plot_config, ax_index=0):
        """
        Plots or updates the thermodynamic diagram on a specified axes.

        This function sets up the axes properties according to the provided plot configuration and
        plots the phase diagram for the specified fluid. It then iterates over all the components in
        the cycle data. For heat exchangers, it plots or updates the processes for both the hot and cold sides.
        For other types of components, it plots or updates the process based on the component's data.

        The function also adjusts the plot limits if it's updating an existing plot,
        ensuring that the axes are scaled correctly to fit the new data.

        Parameters:
        ----------
        ax : matplotlib.axes.Axes
            The axes on which to plot or update the thermodynamic diagram.
        plot_config : dict
            A dictionary containing the plot settings, including variables to be plotted on the x and y axes, and scaling information.
        ax_index : int
            The index of the axes in the figure, used for identifying and updating existing plots.
        """
        # Set up axes properties
        ax.set_xlabel(LABEL_MAPPING[plot_config["x_prop"]])
        ax.set_ylabel(LABEL_MAPPING[plot_config["y_prop"]])
        ax.set_xscale(plot_config["x_scale"])
        ax.set_yscale(plot_config["y_scale"])

        # Plot phase diagram
        cpx.plot_phase_diagram(self.fluid, axes=ax, **plot_config)

        # Plot thermodynamic processes
        for name, component in cycle_data["components"].items():
            # Handle heat exchanger components in a special way
            if component["type"] == "heat_exchanger":
                for side in ["hot_side", "cold_side"]:
                    self._plot_cycle_process(
                        name + "_" + side, plot_config, ax, ax_index=ax_index
                    )
            else:
                self._plot_cycle_process(name, plot_config, ax, ax_index=ax_index)

        # Adjust plot limits if updating
        ax.relim(visible_only=True)
        ax.autoscale_view()

    def _plot_cycle_process(self, name, plot_settings, ax, ax_index=None):
        """
        Creates or updates the plot elements for a specific cycle process on a given axes.

        This method checks if the plot elements for the specified cycle process already exist on the given axes.
        If they exist and an axis index is provided, it updates these elements with new data. Otherwise, it creates
        new plot elements (lines and points) and stores them for future updates.

        Parameters:
        ----------
        name : str
            The name of the cycle process to plot or update.
        plot_settings : dict
            The plot settings dictionary containing settings such as x and y variables.
        ax : matplotlib.axes.Axes
            The axes on which to plot or update the cycle process.
        ax_index : int, optional
            The index of the axes in the figure, used for updating existing plots. If None, new plots are created.
        """

        # Retrieve component data
        x_data, y_data, plot_params = self._get_process_data(
            name, plot_settings["x_prop"], plot_settings["y_prop"]
        )

        # Initialize the dictionary for this axis index if it does not exist
        if ax_index is not None:
            if ax_index not in self.graphics["process_lines"]:
                self.graphics["process_lines"][ax_index] = {}
            if ax_index not in self.graphics["state_points"]:
                self.graphics["state_points"][ax_index] = {}

        # Handle existing plot elements
        if ax_index is not None and name in self.graphics["process_lines"][ax_index]:
            if x_data is None or y_data is None:
                # Hide existing plot elements if data is None
                self.graphics["process_lines"][ax_index][name].set_visible(False)
                self.graphics["state_points"][ax_index][name].set_visible(False)
            else:
                # Update existing plot elements with new data
                self.graphics["process_lines"][ax_index][name].set_data(x_data, y_data)
                self.graphics["state_points"][ax_index][name].set_data(
                    [x_data[0], x_data[-1]], [y_data[0], y_data[-1]]
                )
                self.graphics["process_lines"][ax_index][name].set_visible(True)
                self.graphics["state_points"][ax_index][name].set_visible(True)
        elif x_data is not None and y_data is not None:

            # Prepare kwargs for line and point separately
            line_kwargs = {
                "linestyle": plot_params["linestyle"],
                "linewidth": plot_params["linewidth"],
                "color": plot_params["color"],
                "marker": "none",
                "label": name,
                "zorder": 1,
            }

            point_kwargs = {
                "linestyle": "none",
                "marker": plot_params["marker"],
                "markersize": plot_params["markersize"],
                "markeredgewidth": plot_params["markeredgewidth"],
                "markerfacecolor": plot_params["markerfacecolor"],
                "color": plot_params["color"],
                "zorder": 2,
            }

            # Create new plot elements if data is not None
            x_ends, y_ends = [x_data[0], x_data[-1]], [y_data[0], y_data[-1]]
            (line,) = ax.plot(x_data, y_data, **line_kwargs)
            (points,) = ax.plot(x_ends, y_ends, **point_kwargs)

            # Store the new plot elements
            if ax_index is not None:
                self.graphics["process_lines"][ax_index][name] = line
                self.graphics["state_points"][ax_index][name] = points

    def _get_process_data(self, name, prop_x, prop_y):
        """
        Retrieve thermodynamic data for a specified process of the cycle.

        Parameters:
        ----------
        name : str
            Name of the cycle process.
        prop_x : str
            Property name to plot on the x-axis.
        prop_y : str
            Property name to plot on the y-axis.

        Returns:
        -------
        tuple of (np.ndarray, np.ndarray, str)
            x_data: Array of data points for the x-axis.
            y_data: Array of data points for the y-axis.
            color: Color code for the plot.
        """

        # Normalize public aliases at this plotting boundary. Auxiliary streams
        # must behave identically with canonical names and familiar axis names.
        prop_x = cpx.ALIAS_TO_CANONICAL.get(prop_x, prop_x)
        prop_y = cpx.ALIAS_TO_CANONICAL.get(prop_y, prop_y)

        # Get heat exchanger side data
        if "_hot_side" in name or "_cold_side" in name:
            # Extract the component and side from the name
            if "_hot_side" in name:
                component_name, side_1 = name.replace("_hot_side", ""), "hot_side"
                side_2 = "cold_side"
            else:
                component_name, side_1 = name.replace("_cold_side", ""), "cold_side"
                side_2 = "hot_side"

            data = self.cycle_data["components"][component_name][side_1]
            data_other_side = self.cycle_data["components"][component_name][side_2]
            is_heat_exchanger = True

        else:  # Handle non-heat exchanger components
            data = self.cycle_data["components"][name]
            is_heat_exchanger = False

        # Retrieve data
        is_working_fluid = data["states"]["identifier"] == "working_fluid"
        if (
            not is_working_fluid
            and is_heat_exchanger
            and prop_y == "temperature"
            and prop_x in ["enthalpy", "entropy"]
        ):
            # Special case for heat exchangers in h-T or s-T diagrams
            x_data = data_other_side["states"][prop_x]
            y_data = data["states"][prop_y]
        elif is_heat_exchanger and prop_y == "temperature" and prop_x == "heat_flow":
            # Special case for pinch point diagram
            x_data = data["heat_flow"]
            y_data = data["states"][prop_y]
        elif is_working_fluid:
            # Baseline case for components of the cycle
            x_data = data["states"][prop_x]
            y_data = data["states"][prop_y]

        else:
            # Other cases
            x_data = None
            y_data = None

        # Define plotting specifications
        default_params = {
            "color": "k",  # black
            "linestyle": "-",  # solid line
            "linewidth": 1.25,
            "marker": "o",
            "markersize": 4.5,
            "markeredgewidth": 1.25,
            "markerfacecolor": "w",  # white-filled marker
        }
        plot_params = {**default_params, **data.get("plot_params", {})}

        return x_data, y_data, plot_params

    def _plot_pinch_point_diagram(self, ax, ax_index):
        """
        Plots or updates the pinch point diagram for the thermodynamic cycle's heat exchangers.

        This method visualizes the temperature vs. heat flow rate for each heat exchanger in the cycle.
        The function is capable of plotting the diagram for the first time or updating it with the latest data
        if called subsequently. It uses a sorted approach, beginning with the heat exchanger that has the minimum
        temperature on the cold side and proceeding in ascending order of temperature.

        Parameters:
        ----------
        ax : matplotlib.axes.Axes
            The axes on which to plot or update the pinch point diagram.
        ax_index : int
            The index of the axes in the figure, used to identify and access the specific axes
            for updating the existing plot elements stored in the 'graphics' attribute.

        Notes:
        -----
        The method relies on the 'graphics' attribute of the class to store and update plot elements.
        It handles the creation of new plot elements (lines, endpoints, vertical lines) when first called
        and updates these elements with new data from 'cycle_data' during subsequent calls.
        """

        # Initialize the graphic object dict for this axis index if it does not exist
        if ax_index not in self.graphics["pinch_point_lines"]:
            self.graphics["pinch_point_lines"][ax_index] = {}

        # Set the axes labels
        ax.set_xlabel(LABEL_MAPPING["heat"])
        ax.set_ylabel(LABEL_MAPPING["temperature"])

        # Extract heat exchanger names and their minimum cold side temperatures
        heat_exchangers = [
            (name, min(component["cold_side"]["states"]["temperature"]))
            for name, component in self.cycle_data["components"].items()
            if component["type"] == "heat_exchanger"
        ]

        # Sort heat exchangers by minimum temperature on the cold side
        sorted_heat_exchangers = sorted(heat_exchangers, key=lambda x: x[1])

        # Loop over all heat exchangers
        Q0 = 0.00
        for HX_name, _ in sorted_heat_exchangers:

            # Hot side
            Q_hot, T_hot, plot_params_hot = self._get_process_data(
                HX_name + "_hot_side", "heat_flow", "temperature"
            )
            Q_cold, T_cold, plot_params_cold = self._get_process_data(
                HX_name + "_cold_side", "heat_flow", "temperature"
            )
            # Both temperature arrays are ordered from the cold inlet. Convert
            # hot-side heat released since its inlet to that same coordinate;
            # reversing the values only works for a uniform, symmetric grid.
            Q_hot = Q_hot[0] - Q_hot

            component = self.cycle_data["components"][HX_name]
            # params_hot = component["hot_side"]["plot_params"]
            # props_hot = component["hot_side"]["states"]
            # mass_flow_hot = component["hot_side"]["mass_flow"]
            # Q_hot = (props_hot["h"] - props_hot["h"][0]) * mass_flow_hot
            # T_hot = props_hot["T"]
            # print("old", Q_hot, T_hot)
            # z = zz

            # # Cold side
            # plot_params_cold = component["cold_side"]["plot_params"]
            # props_cold = component["cold_side"]["states"]
            # mass_flow_cold = component["cold_side"]["mass_flow"]
            # Q_cold = (props_cold["h"] - props_cold["h"][0]) * mass_flow_cold
            # T_cold = props_cold["T"]

            # Check if the plot elements for this component already exist
            if HX_name in self.graphics["pinch_point_lines"][ax_index]:
                # Update existing plot elements
                plot_elements = self.graphics["pinch_point_lines"][ax_index][HX_name]
                plot_elements["hot_line"].set_data(Q0 + Q_hot, T_hot)
                plot_elements["cold_line"].set_data(Q0 + Q_cold, T_cold)

                # Update endpoints
                plot_elements["hot_start"].set_data([Q0 + Q_hot[0]], [T_hot[0]])
                plot_elements["hot_end"].set_data([Q0 + Q_hot[-1]], [T_hot[-1]])
                plot_elements["cold_start"].set_data([Q0 + Q_cold[0]], [T_cold[0]])
                plot_elements["cold_end"].set_data([Q0 + Q_cold[-1]], [T_cold[-1]])

                # Update vertical lines
                plot_elements["start_line"].set_xdata([Q0, Q0])
                plot_elements["end_line"].set_xdata([Q0 + Q_hot[-1], Q0 + Q_hot[-1]])
            else:

                # Prepare kwargs for line and point separately
                line_params_hot = {
                    "linestyle": plot_params_hot["linestyle"],
                    "linewidth": plot_params_hot["linewidth"],
                    "color": plot_params_hot["color"],
                    "marker": "none",
                    "zorder": 1,
                }

                point_params_hot = {
                    "linestyle": "none",
                    "marker": plot_params_hot["marker"],
                    "markersize": plot_params_hot["markersize"],
                    "markeredgewidth": plot_params_hot["markeredgewidth"],
                    "markerfacecolor": plot_params_hot["markerfacecolor"],
                    "color": plot_params_hot["color"],
                    "zorder": 2,
                }

                line_params_cold = {
                    "linestyle": plot_params_cold["linestyle"],
                    "linewidth": plot_params_cold["linewidth"],
                    "color": plot_params_cold["color"],
                    "marker": "none",
                    "zorder": 1,
                }

                point_params_cold = {
                    "linestyle": "none",
                    "marker": plot_params_cold["marker"],
                    "markersize": plot_params_cold["markersize"],
                    "markeredgewidth": plot_params_cold["markeredgewidth"],
                    "markerfacecolor": plot_params_cold["markerfacecolor"],
                    "color": plot_params_cold["color"],
                    "zorder": 2,
                }

                # Create new plot elements
                (hot_line,) = ax.plot(Q0 + Q_hot, T_hot, **line_params_hot)
                (cold_line,) = ax.plot(Q0 + Q_cold, T_cold, **line_params_cold)

                # Create endpoints
                (hot_1,) = ax.plot(Q0 + Q_hot[0], T_hot[0], **point_params_hot)
                (hot_2,) = ax.plot(Q0 + Q_hot[-1], T_hot[-1], **point_params_hot)
                (cold_1,) = ax.plot(Q0 + Q_cold[0], T_cold[0], **point_params_cold)
                (cold_2,) = ax.plot(Q0 + Q_cold[-1], T_cold[-1], **point_params_cold)

                # Create vertical lines
                param = {"color": "black", "linestyle": "-", "linewidth": 0.75}
                start_line = ax.axvline(x=Q0, zorder=1, **param)
                end_line = ax.axvline(x=Q0 + Q_hot[-1], zorder=1, **param)

                # Store new plot elements
                self.graphics["pinch_point_lines"][ax_index][HX_name] = {
                    "hot_line": hot_line,
                    "cold_line": cold_line,
                    "hot_start": hot_1,
                    "hot_end": hot_2,
                    "cold_start": cold_1,
                    "cold_end": cold_2,
                    "start_line": start_line,
                    "end_line": end_line,
                }

            # Update abscissa for the next heat exchanger
            Q0 += Q_hot[-1]

        ax.set_xlim(left=0 - Q0 / 50, right=Q0 + Q0 / 50)

    def _get_diagram_default_settings(self, plot_config):
        """
        Merges user-provided plot settings with default settings.

        Parameters
        ----------
        user_settings : dict
            A dictionary of user-defined plot settings.

        Returns
        -------
        dict
            A dictionary containing the merged plot settings.
        """
        default_settings = {
            "x_prop": "entropy",
            "y_prop": "temperature",
            "x_scale": "linear",
            "y_scale": "linear",
            "plot_saturation_line": True,
            "plot_critical_point": True,
            "plot_quality_isolines": False,
            "plot_pseudocritical_line": False,
            "plot_triple_point_vapor": False,
            "plot_triple_point_liquid": False,
            "plot_spinodal_line": False,
            # Add other default settings here
        }

        # Combine the global fluid settings with the "plots" settings
        fluid_settings = self.plot_settings.get("fluid", {})
        fluid_settings = {} if fluid_settings is None else fluid_settings
        plot_config = plot_config | fluid_settings

        # Merge with default values
        settings = default_settings | plot_config
        for axis in ("x_prop", "y_prop"):
            settings[axis] = cpx.ALIAS_TO_CANONICAL.get(settings[axis], settings[axis])
        return settings

    def save_data_to_excel(self, filename="performance.xlsx"):
        """
        Exports the cycle performance data to Excel file
        """

        if "cycle_data" not in self.__dict__:
            raise RuntimeError("Evaluate fitness after configuring or updating the problem before exporting cycle results.")

        # Define variable map
        variable_map = {
            "fluid_name": {"name": "fluid_name", "unit": "-"},
            "temperature": {"name": "temperature", "unit": "K"},
            "pressure": {"name": "pressure", "unit": "Pa"},
            "density": {"name": "density", "unit": "kg/m3"},
            "quality_mass": {"name": "quality", "unit": "-"},
            "compressibility_factor": {"name": "compressibility_factor", "unit": "-"},
            "internal_energy": {"name": "internal_energy", "unit": "J/kg"},
            "enthalpy": {"name": "enthalpy", "unit": "J/kg"},
            "entropy": {"name": "entropy", "unit": "J/kg/K"},
            "isobaric_heat_capacity": {"name": "isobaric_heat_capacity", "unit": "J/kg/K"},
            "isochoric_heat_capacity": {"name": "isochoric_heat_capacity", "unit": "J/kg/K"},
            "heat_capacity_ratio": {"name": "heat_capacity_ratio", "unit": "-"},
            "speed_of_sound": {"name": "speed_of_sound", "unit": "m/s"},
            "viscosity": {"name": "dynamic_viscosity", "unit": "Pa*s"},
            "conductivity": {"name": "thermal_conductivity", "unit": "W/m/K"},
            "superheating": {"name": "superheating_degree", "unit": "K"},
            "subcooling": {"name": "subcooling_degree", "unit": "K"},
        }

        # Initialize a list to hold all rows of the DataFrame
        data_rows = []

        # Prepare the headers and units rows
        headers = ["state"]
        units_row = ["units"]

        for key in variable_map:
            headers.append(variable_map[key]["name"])
            units_row.append(variable_map[key]["unit"])

        def export_value(state, key):
            # Optional diagnostics can be disabled or unsupported by a model.
            # Export a blank numeric cell rather than inventing a zero or
            # performing extra property flashes while writing the workbook.
            try:
                # Excel writers do not recognize zero-dimensional JAX arrays
                # as numbers. Materialize scalar cells only at this I/O edge.
                return np.asarray(state[key]).item()
            except (cpx.UnsupportedPropertyError, cpx.PropertyUnavailableError):
                return np.nan

        # Iterate over each component in the dictionary
        for component_name, component in self.cycle_data["components"].items():
            if component["type"] == "heat_exchanger":
                # Handle heat exchanger sides separately
                for side in ["hot_side", "cold_side"]:
                    # Append the data for state_in and state_out to the rows list
                    state_in = component[side]["state_in"]
                    state_out = component[side]["state_out"]
                    data_rows.append(
                        [f"{component_name}_{side}_in"]
                        + [export_value(state_in, key) for key in variable_map]
                    )
                    data_rows.append(
                        [f"{component_name}_{side}_out"]
                        + [export_value(state_out, key) for key in variable_map]
                    )
            else:
                # Handle non-heat exchanger components
                # Append the data for state_in and state_out to the rows list
                state_in = component["state_in"]
                state_out = component["state_out"]
                data_rows.append(
                    [f"{component_name}_in"]
                    + [export_value(state_in, key) for key in variable_map]
                )
                data_rows.append(
                    [f"{component_name}_out"]
                    + [export_value(state_out, key) for key in variable_map]
                )

        # Create a DataFrame with data rows
        df = pd.DataFrame(data_rows, columns=headers)

        # Insert the units row
        df.loc[-1] = units_row  # Adding a row
        df.index = df.index + 1  # Shifting index
        df = df.sort_index()  # Sorting by index

        # Prepare energy_analysis data
        df_2 = pd.DataFrame(
            list(self.cycle_data["energy_analysis"].items()),
            columns=["Parameter", "Value"],
        )

        # Export to Excel
        # filename = os.path.join(self.out_dir, filename)
        with pd.ExcelWriter(filename, engine="openpyxl") as writer:
            df.to_excel(writer, index=False, sheet_name="cycle_states")
            df_2.to_excel(writer, index=False, sheet_name="energy_analysis")
