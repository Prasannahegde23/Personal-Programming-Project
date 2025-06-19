Readme file for Bond failure analyses test:

This folder contains necessary files and folders to run a Bond failure analyses test. 
- with_hydrogen_and_crack and without_hydrogen_no_crack folers contain files ready for respective testing.
- with_hydrogen.txt and without_hydrogen.txt files are the results of the test.
- Bond_failure_analyses.py python script is used to generate and asses the respective bond damage plot.

Use the following command in your terminal to run the python script:
python Bond_failure_analyses_plot.py

To run Peridigm:
Peridigm <Input_file.yaml>

Note:
- Since the Output.txt file opens in append mode, 
it needs to be deleted each and every time for the new run.
- Turn off the boundary condition test in .yaml file under material section