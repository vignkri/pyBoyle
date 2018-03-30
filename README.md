# Boyle - Python

Computation agent for Boyle made using python-scipy. The model is built based on the Anaerobic Digestion model for computing the gaseous fraction of the methane produced from the digestion process.

## Version Changelog

**Current Version**: 0.4

The feature list of this version are:

1. Configuration file using YAML to set up the program without complicated commands and interactivity.
2. Simulation manager to run the simulation within a secure zone.
3. Exporter is updated to output a hdf5 file. The output data structure is given in the section [#Notes on Data]
4. The logger now handles independent module based messages for debugging.
5. IOTools has been re-factored to perform both input and output settings.
6. Boyle can now be called similar to a CLI application with a path function for specific simulation models. 

## How to Use

Currently, this is not a package nor an application. This is just a collection of scripts designed to run a simulation model. In order to run the model, follow the steps:

1. Make sure Python 3.x is installed in your computer.
2. Clone the repository using `git clone <url>`.
3. Navigate to the root folder and run `pipenv install` to set up the python virtual environment or the tool of your choice to setup your virtual environment.
4. Activate the environment by `pipenv shell`
6. Navigate to the boyle folder and run `python boyle.py -path "./simulation.yaml"` to run the simulation.
    1. Make sure the `simulation.yaml` file points to the correct data source for the simulation to run properly.

## Notes on Data

1. The input data is in `grams/liter` thereby providing the density.
2. The output data structure:
    1. Contains three different groups: *metadata, debug, solution, headers*
    2. `Metadata` contains specific information on when the experiment / simulation was performed, the name of the experiment, and the time it finished.
    3. `Headers` contains the header items for parsing the outputs properly.
    4. `Debug` contains the output of results without removing are adjusting specific output conditions. This allows for debugging resultant data.
    5. `Solution` contains the final solution data where some columns are curtailed so as to have an untouched output file.
