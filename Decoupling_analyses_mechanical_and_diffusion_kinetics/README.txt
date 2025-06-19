Readme file for Decoupling analyses for mechanical and diffusion kinetics:

This folder contains necessary files and folders to run a 
Decoupling analyses for mechanical and diffusion kinetics test: 
- "Hydrogen_only"" and "Mechanical_only" folers contain files ready for respective testing.
- Hydrogen_only.txt and Mechanical_only.txt files are the results of the test placed under Analysis.
- Post_processing_PDHE.py python script is used to generate and asses
 final contour plot of respective cases. These contour plots obtained, can be used for 
 assessing Damage, Displacements and Hydrogen concentration.

To run the python script, use the following command in your terminal:
python Post_processing_PDHE.py <file_name>.txt

To run Peridigm:
Peridigm <Input_file.yaml>

Note:
- Since the Output.txt file opens in append mode, 
it needs to be deleted each and every time for the new run.
- Turn off the boundary condition test in .yaml file under material section