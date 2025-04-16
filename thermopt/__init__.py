# # Highlight exception messages
# # https://stackoverflow.com/questions/25109105/how-to-colorize-the-output-of-python-errors-in-the-gnome-terminal/52797444#52797444
# try:
#     import IPython.core.ultratb
# except ImportError:
#     # No IPython. Use default exception printing.
#     pass
# else:
#     import sys
#     sys.excepthook = IPython.core.ultratb.FormattedTB(color_scheme='linux', call_pdb=False)

# Import submodules
from .math import *
# from .plot_functions import *
from .config import *

# Import subpackages
from pysolver_view import *
from .properties import *
from .utilities import *


from .cycles import *
from .properties import *

# Set plot options
set_plot_options()


# Package info
__version__ = "0.1.0"
PACKAGE_NAME = "thermopt"
URL_GITHUB = "https://github.com/turbo-sim/thermopt"
URL_DOCS = "https://turbo-sim.github.io/thermopt/"
URL_DTU = "https://thermalpower.dtu.dk/"
BREAKLINE = 80 * "-"


def print_banner():
    """Prints a banner."""
    banner = r"""
      ________                        ____        __ 
     /_  __/ /_  ___  _________ ___  / __ \____  / /_
      / / / __ \/ _ \/ ___/ __ `__ \/ / / / __ \/ __/
     / / / / / /  __/ /  / / / / / / /_/ / /_/ / /_  
    /_/ /_/ /_/\___/_/  /_/ /_/ /_/\____/ .___/\__/  
                                       /_/           
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


