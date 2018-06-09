__version__ = "0.6.0"

# Set up imports

from boyle.core.generic import Dataset, SimulationResult
from boyle.core import save, load
from boyle.manager import Manager

# Set up imports from api in submodules
from boyle.tools.api import *
from boyle import preprocessing
