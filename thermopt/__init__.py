# from .cycles.components import *

# Import submodules
from .math import *
# from .plot_functions import *
from .config_validation import *

# Import subpackages
from .pysolver_view import *
from .properties import *
from .utilities import *


from .cycles import *
from .properties import *

# from .brayton_recuperated import *
# from .brayton_split_compression import *



# Set plot options
set_plot_options()


# Package info
__version__ = "0.0.1"
PACKAGE_NAME = "thermopt"
URL_GITHUB = "https://github.com/turbo-sim/thermopt"
URL_DOCS = "https://turbo-sim.github.io/thermopt/"
URL_DTU = "https://thermalpower.dtu.dk/"
BREAKLINE = 80 * "-"


def print_banner():
    """Prints a banner."""
    banner = """
        ______ __  __ ____   ____   ____   ______ __    ____  _       __
       /_  __// / / // __ \ / __ ) / __ \ / ____// /   / __ \| |     / /
        / /  / / / // /_/ // __  |/ / / // /_   / /   / / / /| | /| / / 
       / /  / /_/ // _, _// /_/ // /_/ // __/  / /___/ /_/ / | |/ |/ /   
      /_/   \____//_/ |_|/_____/ \____//_/    /_____/\____/  |__/|__/  
    """
    print(BREAKLINE)
    print(banner)
    print(BREAKLINE)
    # https://manytools.org/hacker-tools/ascii-banner/
    # Style: Sland


def print_package_info():
    """Prints package information with predefined values."""

    info = f""" Version:       {__version__}
 Repository:    {URL_GITHUB}
 Documentation: {URL_DOCS}"""
    print_banner()
    print(BREAKLINE)
    print(info)
    print(BREAKLINE)


