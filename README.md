# SA_BIM2SIM
Influence of parameters on building simulations of Building Information Modelling (BIM) models using Reduced Order Model (ROM) approach

Requirements: All files in SA_BIM2SIM, Dymola-Python interface, ebcpy_modelica_simres.py and required packages in scrpits, exported data in form of pickle file of BIM model using BIM2SIM toolchain 

Run automatic.py
Enter path where pickle file is kept in computer
Enter name of pickle file to be selected for analysis
Enter path where models need to be exported for simulations ( 180 models for Morris analysis and 5000 models for Sobol analysis)

Wait till simulation completes

Select multiprocessing of simulations or single model simulations
If chosen multiprocessing, enter number of simulations to run simultaneously

Note: Change directory in scripts as per your computer for libraries locations (like Aixlib, BIM2SIM)
