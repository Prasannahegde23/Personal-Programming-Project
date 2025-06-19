Readme file for Mesh convergence test:

This folder contains necessary files to run a Mesh convergence test. 

Use the following command in your terminal:
Peridigm <Input_file.yaml>

Displacements_output and Hydrogen_output folders contain the output files for 
all considered mesh densities and respective python scripts are used to plot and asess the 
mesh convergence.

To run the python script:
python <plotting_file.py>

After running the Peridigm input file, you will see the results under "Output.txt" file.

Note:
- Since the Output.txt file opens in append mode, 
it needs to be deleted each and every time for the new run.
- Turn off the boundary condition test in .yaml file under material section