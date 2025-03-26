# ThermoOpt

``ThermoOpt`` is a Python package for the modeling and optimization of thermodynamic cycles.

📚 **Documentation**: [https://turbo-sim.github.io/thermopt/](https://turbo-sim.github.io/thermopt/) *(under construction)*  
📦 **PyPI package**: [https://pypi.org/project/thermopt/](https://pypi.org/project/thermopt/)


## 🚀 User installation (via PyPI)

If you just want to use ``ThermoOpt``, the easiest way is to install it from PyPI:

```bash
pip install thermopt
```


You can then verify the installation with:

```bash
python -c "import thermopt; thermopt.print_package_info()"
```


## 🛠️ Developer installation (from source with Poetry)

This guide walks you through installation for development using `Poetry`, which manages both dependencies and virtual environments automatically.

1. **Install Poetry**  
   Follow the official guide: [Poetry Installation](https://python-poetry.org/docs/#installation)  
   Then verify the installation:
   ```bash
   poetry --version
   ```

2. **Clone the repository**
   ```bash
   git clone https://github.com/turbo-sim/thermopt
   cd thermopt
   ```

3. **Install dependencies and create a virtual environment**  
   This installs all required packages listed in `pyproject.toml`:
   ```bash
   poetry install
   ```

4. **Activate the virtual environment**
   ```bash
   poetry shell
   ```
   To deactivate:
   ```bash
   exit
   ```

5. **Verify the installation**  
   Run the following inside the Poetry shell or with `poetry run`:
   ```bash
   poetry run python -c "import thermopt; thermopt.print_package_info()"
   ```

6. **Install additional packages**  
   To add a runtime dependency:
   ```bash
   poetry add <package-name>
   ```
   To add a development-only dependency:
   ```bash
   poetry add --dev <package-name>
   ```


## 📂 Examples

The [examples](examples) directory contains a variety of ready-to-run thermodynamic cycle cases, covering different working fluids and applications.

Each example:
- Is defined in a `.yaml` input file
- Is executed via a corresponding `run_optimization.py` script
- Outputs results in a subdirectory called `results/`

To run any example, navigate to the corresponding subfolder and execute the optimization script.