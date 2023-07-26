# Dynamics of Vehicles - HW1 -Tyre fitting
The current repository contains the code developed during Master's Degree in Mechatronics Engineering, for the project of the course "Dynamics of Vehicles" (Biral Francesco) 

## Assignment 1
Tyre fitting using MF96 formula ([_Pacejka_](https://en.wikipedia.org/wiki/Hans_B._Pacejka) Magic Formula 96).
In the first assignment each team had to fit the MF96 tyre model coefficients on a dataset of measured tyre forces for a F-SAE tyre. The specific tyre used for the fitting was "Hoosier 18.0 x 6.0 - 10 R25B". The main purpose was to obtain acceptable fitting results for different cases: pure longitudinal/lateral forces, combined longitudinal/lateral forces and self aligning moment.

### 
In the folder MATLab we can find the main code developed during April 2023 to finish the work relative to this first assignment.
In particular there is the [main file](MATLab/tyre_lib/main_tyre_data_analysis.m) contating lots of code recalling lot of files in tyre_lib to deal with the different cases for fitting the forces and moments.
