Readme file for Installation files:

This folder contains necessary shell scripts for directory setup 
and peridigm and meshio installation into your local node. Detailed explanation of 
how to install peridigm using these files are given in the report.

For Peridigm, first use Directory_setup.sh and then make-peridigm.sh
Use the following command in your terminal to change permissions and run the scripts:
chmod u+x <bash_script.sh>
./<bash_script.sh>

Place the content present in "Content_to_be_placed_bashrc_zshrc.txt" in your 
bashrc or zshrc file in order to load dependecies.

After successfull installation check:
which Peridigm
Peridigm --version
Peridigm <Input_file.yaml>
