Readme file for Unit test:

This folder contains necessary files to run a Unit test.

Use the following command in your terminal:
g++ -std=c++17 -O2 Unit_test_PDHE.cpp PDHE_element_routine.cxx PDHE_material_routine.cxx PDHE_boundary_effects.cxx -o Unit_test_PDHE -lstdc++fs

After compiling, run the generated shell-script to botain log file
./Unit_test_PDHE