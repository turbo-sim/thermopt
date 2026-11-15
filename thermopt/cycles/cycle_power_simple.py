from . import cycle_power_recuperated

def evaluate_cycle(
    variables,
    parameters,
    constraints,
    objective_function,
    fluids=None,
):
    
    output = cycle_power_recuperated.evaluate_cycle(
        variables,
        parameters,
        constraints,
        objective_function,
        recuperated=False,
        fluids=fluids,
    )

    return output
