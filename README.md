# Boyle - Python

Computation agent for Boyle made using python-scipy. The model is built based on the Anaerobic Digestion model for computing the gaseous fraction of the methane produced from the digestion process.

## Version Changelog

**Current Version**: 0.6.2

The current version of the simulation engine is designed as a package that can be imported to any IDE or Jupyer-Notebook in order to facilitate flexibility in operation of the simulations. This allows for flexible assignment of outputs and settings in order to perform the simulation. The documentation is still under review until the documentations are up, a series of sample scripts and tools will be updated for general use cases.

## Notes on Data

1. The input data is in `grams/liter` thereby providing the density.
2. The output data structure:
    1. Contains three different groups: *metadata, debug, solution, headers*
    2. `Metadata` contains specific information on when the experiment / simulation was performed, the name of the experiment, and the time it finished.
    3. `Headers` contains the header items for parsing the outputs properly.
    4. `Debug` contains the output of results without removing are adjusting specific output conditions. This allows for debugging resultant data.
    5. `Solution` contains the final solution data where some columns are curtailed so as to have an untouched output file.
