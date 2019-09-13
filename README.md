# ratiometric-gpcr-particle-sims
Particle-based simulation and analysis code for "Ratiometric GPCR signaling enables directional sensing in yeast".

The simulations were performed on Longleaf, a UNC compute cluster, using [Smoldyn](www.smoldyn.org) v2.56 and analyzed in [MATLAB](https://www.mathworks.com/products/matlab.html) 2017b. Representative scripts to generate, process, and analyze main text figures are provided. Additionally, .csv files containing data from Figures 5g, 6b, and 8e are provided in Representative_Analyzed_Data. A more complete set of .csv files for all the particle-based simulation figure panels in the paper is provided in __Data_For_Figures__.

# Workflow
The folders and files provided are meant to reflect the typical workflow for analyzing the data, and are not all intended to be usable out-of-the-box. The following assumes Smoldyn and MATLAB are set up appropriately. Simulations can be run on a desktop computer, with run times on the order of day(s); the UNC Longleaf compute cluster was used to perform simulations in ["embarrassingly parallel"](https://en.wikipedia.org/wiki/Embarrassingly_parallel) fashion.

1. Generate configuration files for Smoldyn, specifying parameters of interest.
2. Run the configuration files using Smoldyn.
3. Analyze the output coordinate data files in MATLAB.

### Config_File_Generation
The function `generate_screen()` sets parameters to simulate. `generate_screen()` calls `create_cfg()` to produce Smoldyn configuration files for each parameter set of interest. The resulting configuration files are plain text files, labeled as `.cfg` files. `create_cfg()` also generates a Bash file to launch Longleaf jobs for simulating the entire set of `.cfg` files created in a single call to `generate_screen()`.

### Representative_Config_Files
Representative `.cfg` files. A single `.cfg` file is reasonable to run on a local machine by entering `smoldyn myconfig.cfg` at the command line; (see the [Smoldyn](www.smoldyn.org) webpage for installation and use details); typical simulations took on the ~1 day on a Mac with 3.4 GHz Intel processor.

### Representative_Raw_Data
The output from running Smoldyn with a `.cfg` file is an xyz coordinate data file, organized by particle type and timestep, denoted `.xyz`. Representative files `.xyz` files are provided in this folder.

### Example_Analysis_Scripts
Example scripts for analyzing `.xyz` files. The main quantity of interest is the angle between the active G protein center of mass and the gradient.

### Representative_Analyzed_Data
CSV files containing analyzed data corresponding to Figures 5g, 6b, and 8e are provided in Representative_Analyzed_Data. Each file contains 50 realizations of the labeled condition, over 600 seconds of simulation.
