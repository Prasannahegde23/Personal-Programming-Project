#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 08:14:03 2025

@author: prasannahegde
"""

import meshio
import sys

# 1) Read the Gmsh file
input_msh_file = sys.argv[1]
output_exo_file = sys.argv[2]
mesh = meshio.read(input_msh_file)

# 2) Create and open a text file for writing
with open(output_exo_file, "w") as f:
    # Optional: You can write a header if you like, but Peridigm's TextFile discretization
    # typically does not need one. If you do include it, remove it or comment it out:
    # f.write("x y z Volume BlockID\n")

    # 3) For each node, write out x,y,z plus placeholders for Volume & BlockID
    for i, point in enumerate(mesh.points):
        x, y, z = point*1e-3
        volume = 1.420393942*1e-6  # placeholder or computed
        blockID = 1      # placeholder or computed
        f.write(f"{x} {y} 0.0 {volume} {blockID}\n")
