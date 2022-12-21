# SMIS
## Single Molecule Imaging Simulator

SMIS is a Matlab-based simulation software that enables to simulate a large variety of single-molecule fluorescence imaging experiments in Widefield mode. 

The novelty of SMIS is the advanced description of the fluorophores used in the simulations, notably taking into account their spectral and photophysical characteristics.

Examples of simulations that can be performed include PALM, dSTORM, sptPALM, PAINT, in 2D or 3D. 

Multicolor experiments can be simulated with unlimited number of colors. 

Complex laser excitation schemes can be simulated with unlimited number of lasers. 

Complex diffusion patterns of the fluorophores can be simulated with unlimited number of diffusion states.

SMIS outputs .tif image stacks and ground truth .mat data.

## =======================================================
Copyright (C) 2022 Dominique Bourgeois
E-mail: dominique.bourgeois@ibs.fr

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
## ======================================================= 

# Reference
Bourgeois, D.; Single Molecule Imaging Simulations with Advanced Fluorophore Photophysics. https://www.biorxiv.org/content/10.1101/2022.06.14.496133v2

## Reported bugs (to be fixed in the next SMIS release)
None for now ...

# Installation
The software can be used either as a standalone application for Windows, MacOS or Linux, or as a MATLAB app. 
SMIS was developed under Windows. Proper running of SMIS under MacOS or Linux has not been thoroughly checked.

## Downloading SMIS

For Windows and Linux, use the default SMIS branch called main_LINUX_WINDOWS; 
Click on the green button “Code”, and select “download ZIP”. 
This will download the file: “SMIS-main_LINUX_WINDOWS.zip”

For Mac OS, switched to the SMIS branch called MACOS; 
Click on the green button “Code”, and select “download ZIP”. 
This will download the file: “SMIS-MACOS.zip”

Finally, unzip the SMIS .zip file in your preferred directory.

### Standalone SMIS 
To run the standalone SMIS, the Matlab runtime must be installed on your computer. The runtime with proper version should be installed. To download the Matlab runtime, go to https://fr.mathworks.com/products/compiler/matlab-runtime.html

=> Windows:
Make sure Matlab runtime 2022a (9.12) is installed on your computer.
The SMIS executable is found in: SMIS/STANDALONE/DISTRIBUTE/WINDOWS
To execute SMIS, double-click on “SMIS.exe”

=> Linux:
Make sure Matlab runtime 2021b (9.11) is installed on your computer.
The SMIS executable is found in: SMIS/STANDALONE/DISTRIBUTE/LINUX
To execute SMIS, open a terminal, move to the directory where SMIS is installed and type at the prompt: “./run_SMIS.sh <mcr_directory>”, where <mcr_directory> is the location of the Matlab runtime. 

=> Mac OS:
Make sure Matlab runtime 2022b (9.13) is installed on your computer.
The SMIS executable is found in: SMIS/STANDALONE/DISTRIBUTE/MACOS
To execute SMIS, open a terminal, move to the directory where SMIS is installed and type at the prompt: “./run_SMIS.sh <mcr_directory>”, where <mcr_directory> is the location of the Matlab runtime. 

### Using the SMIS app with Matlab
To properly run SMIS with Matlab, you need the Image Processing Toolbox, and preferably the Parallel Computing Toolbox (not compulsory).

In Matlab, go to the APPS tab, and click on “Install App”. 

=> For Windows select the “SMIS_Windows.mlappinstall” located in SMIS/APP 

=> For Linux select the “SMIS_Linux.mlappinstall” located in SMIS/APP 

=> For Mac OS select the “SMIS.mlappinstall” located in SMIS/APP 

## SMIS User manual

Please look at the user manual for detailed instructions how to use SMIS.

In the SMIS GUI, interactive help is also available by moving the mouse to the desired field.

### Getting familiar with SMIS: 

To start learning about SMIS, it is advisable to load simulation examples by clicking on: “Load simulations”. Then navigate throughout the different SMIS submenus. Enter a proper output directory for the simulation and click on “Launch simulation”.

### Main steps to design a simulation : 
1/ If there is a similar simulation available, load it with « Load simulation parameters ».

2/ Enter the name of your simulation in the « SMIS Simulation Title ».

3/ Choose whether this is a 2D or 3D simulation with « 3D Off-On ».

4/ Set the number of different fluorophores (e.g. for multicolor experiments) you want to use in « Number of fluorophores ».

5/ Load the virtual samples for each defined fluorophore via the « Load virtual sample » menu. You can create virtual samples with the menu “Create virtual samples”.

6/ Set up labeling and photophysics for each fluorophore with the « Setup Fluorophore » menu. In this menu, you can load existing fluorophores, or create your own fluorophore by entering “Define new fluorophore”.

7/ Decide whether this is « sptPALM », « qPALM » or « FRET » experiment. The fluorophores and virtual samples must have been chosen accordingly.

8/ Defined the number of lasers to be used and set up each laser with the « Setup Laser » menu.

9/ Define the final frame size of the detector (and output stack) by setting the « Binning factor »

10/ Define the « Number of frames », « Raster size », « Frametime » and « Time between frame time ».

11/ Define eventual « Sample drift ».

12/ Decide whether this is a single or two-channel experiment with the « Number of acquisition channels » toggle.

13/ Set up the microscope « Objective and PSF » parameters.

14/ Set up the « Emission Filters », and eventually the parameters of the « EMCCD camera » (such as the EMCCD gain).

15/ Define the « Fluorescence background ».

16/ Choose the « Output Directory » and « Stack File Name ».

17/ Save your simulation with the « Save simulation parameters » button. Saving can be done at any time during the process.

18/ Finally launch the simulation with the « Launch simulation » button (single molecule mode) or « Launch ensemble simulation » (ensemble mode).

