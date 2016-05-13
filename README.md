# 1DMultilayerAdditiveManufacture
The code is written in MATLAB and aims to simulate a multi-layer thermal problem using linear FEM.
The problem we want to solve is a 1D thermal problem with moving right boundary condition. The process is solved 
twice: firstly on a global coarse mesh and secondly on a locally (on the coarse right-end element) refined mesh.
The refined mesh gets the boundary conditions from the coarse mesh and solve the local problem independently.
Two strategy are then employed to exchange the information from the local to the global mesh:

1) L2 projection onto the coarse mesh

2) X-FEM enrihment of the global mesh basis after a POD training process to generate the enrichment modes


Tests:
The code is tested using the NAFEM Benchmark where at x=0.08 the reference solution is 36.60°C after 32sec using 10 elements.
The result with the present code is 35.64°C, returning a relative error of 0.026229508.
