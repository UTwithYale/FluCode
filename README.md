# FluCode - UT Austin / Yale
The model can be run by executing the main file: Jul17_Step2_sims_humidity.mlx

The resulting .mat will be saved in folder: /augMats.

This file creates .mat files. Each file contains a single simulation for each of the 32 scenarios.

To format the results in the common template format the file the script toTemplate_main.m should be executed. It parses the .mat files located in the folder mats and exports csv results in the folder ExcelOuput.
