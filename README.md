# Boyle - Python

Computation agent for Boyle made using python-scipy. The model is built based on the Anaerobic Digestion model for computing the gaseous fraction of the methane produced from the digestion process.

## Version Changelog

**Current Version**: 0.2

The feature list of this version are:

1. Configuration file using YAML to set up the program without complicated commands and interactivity.
2. Simulation manager to run the simulation within a secure zone.
3. Exporter is updated to output a single pickle file instead of CSVs.
4. The logger now handles independent module based messages for debugging.

## How to Use

Currently, this is not a package nor an application. This is just a collection of scripts designed to run a simulation model. In order to run the model, follow the steps:

1. Make sure Python 3.x is installed in your computer.
2. Clone the repository using `git clone <url>`.
3. Navigate to the root folder and run `virtualenv .venv` to set up the python virtual environment or the tool of your choice to setup your virtual environment.
4. Make sure the virtual environment is active.
5. Install the requirements using `pip install -r requirements.txt` from the root folder.
6. Navigate to the boyle folder and run `boyle.py` to run the simulation.
