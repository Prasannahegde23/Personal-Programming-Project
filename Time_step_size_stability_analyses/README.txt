Readme file for Time step-size analyses test:

This folder contains necessary files and folders to run a Time step-size analyses test. 
- Folders starting with name "critic_<alpha_value>" contains files to run the test.
- "stable" folder contains files where time step size is taken according to CFL-stability.
- "Analysis" folder contains all output files and a python script to plot the convergence test.

Use the following command in your terminal to run the python script:
python For_euler_scheme.py

To run Peridigm:
Peridigm <Input_file.yaml>

Note:
- Since the Output.txt file opens in append mode, 
it needs to be deleted each and every time for the new run.
- Turn off the boundary condition test in .yaml file under material section