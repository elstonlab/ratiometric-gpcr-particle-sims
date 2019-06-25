# ratiometric-gpcr-particle-sims
Particle-based simulation and analysis code for "Ratiometric GPCR signaling enables directional sensing in yeast".

The simulations were performed on Longleaf, a UNC compute cluster, using [Smoldyn](www.smoldyn.org) v2.56 and analyzed in [MATLAB](https://www.mathworks.com/products/matlab.html) 2017b. Representative scripts to generate, process, and analyze main text figures are provided. Additionally, .csv files containing data from Figures 5g, 6b, and 8e are provided in Representative_Analyzed_Data.

# Workflow
The folders and files provided are meant to reflect the typical workflow for analyzing the data, and are not all intended to be usable out-of-the-box. The following assumes Smoldyn and MATLAB are set up appropriately. Simulations can be run on a desktop computer, with run times on the order of day(s); the UNC Longleaf compute cluster was used to perform simulations in ["embarrassingly parallel"](https://www.mathworks.com/products/matlab.html) fashion.

## Folder details
### Config_File_Generation
Sets parameters to simulate in `generate_screen()`, which calls `create_cfg()` to produce Smoldyn configuration files. These are plain text configuration files that I label as `.cfg` files. This function also generates a Bash file to launch Longleaf jobs for simulating the entire set of `.cfg` files created in a single call to `generate_screen()`.

### Representativ_Config_Files
Representative `.cfg` files. A single `.cfg` file is reasonable to run on a local machine (see the [Smoldyn](www.smoldyn.org) webpage for installation and use); typical parameters took on the order of 1 day on a Mac with 3.4 GHz Intel processor. 

### Representative_Raw_Data
The output from running Smoldyn with the `.cfg` files is a set of coordinate data organized by particle type and timestep, denoted `.xyz`. Representative files are provided in this folder.

### Example_Analysis_Scripts
Example scripts for analyzing `.xyz` files.

### Representative_Analyzed_Data
CSV files containing data from Figures 5g, 6b, and 8e are provided in Representative_Analyzed_Data. Each file contains 50 realizations of the labeled condition, over 600 seconds of simulation.
