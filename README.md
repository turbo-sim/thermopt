
# ThermoOpt

``ThermoOpt`` is a python package for the modeling and optimization of thermodynamic cycles.


## Quick installation guide (from Git repository)

This guide will walk you through installing ``Thermopt`` using ``Poetry``, which manages dependencies and virtual environments efficiently.

### 1. Install Poetry

If you haven't installed Poetry yet, follow the official guide:  
[Poetry Installation](https://python-poetry.org/docs/#installation)


Verify the installation:
```sh
poetry --version
```

### 2. Clone the Repository

If you haven't cloned the ``ThermOpt`` repository yet, do so:
```sh
git clone https://github.com/turbo-sim/thermopt
cd thermopt
```

### 3. Create a Virtual Environment and Install Dependencies

Poetry will automatically create a virtual environment and install all dependencies specified in `pyproject.toml`:
```sh
poetry install
```

### 4. Activate the Virtual Environment

To activate the virtual environment, run:
```sh
poetry shell
```

To deactivate it, simply type:
```sh
exit
```

### 5. Verify the Installation

Run the following command inside the Poetry shell or using `poetry run`:
```sh
poetry run python -c "import thermopt; thermopt.print_package_info()"
```
If successful, this will display the Thermopt banner and package details.

### 6. Installing Additional Packages

To install new dependencies:
```sh
poetry add <package-name>
```

For development dependencies:
```sh
poetry add --dev <package-name>
```

---
