
# Description #
This code solves the 1D dynamic compaction of a two-phase granular medium using the algorithm described in "Simple and efficient relaxation methods for interfaces separating compressible fluids, cavitating flows and shocks in multiphase mixtures" (Saurel 2018). The test cases from the paper can be run in `runfiles/compare_saurel.m`.

- First set the path of your repository in the variable `maindir` in the file `main/runfiles/pathfile.m`.

For a 1D dynamic compaction, a velocity boundary condition is imposed on the left boundary, and the mesh is dynamically updated. 
* Set initial conditions and variables in `initial_simple.m`.
	- Set `plotit = 1` to plot results during simulation.
	- Set `savefile = 1` to save the results.
	- Set number of monte carlo realizations in `mc_count`.

* For a single realization, run `simple_compact.m`.
* For a Monte Carlo simulation run `mc_compact.m`. Set the variable `mc_count` for the number of realizations. 
* For constant density, check files and folders with suffix `constrho`.

- Data is stored in `results/data/` file. 
- Scripts for reading data, analyzing it and plotting are in `results/script/`. 
- `readtxtfile.m` reads the data of a single realization
- `mc_read.m` with input the name of the folder containing the MC output files to read a whole Monte Carlo with multiple realizations.
- `mc_analysis.m` reads MC data and plots PDF, CDF and temperature ditributions.

- The folder `helpfunc` contains helper functions used in the main code; typically small functions.
- The folder `solver` contains the main functions of the code.
- The function `solver1D.m` contains the main solver.

[//]: <> (


- Fix:
	- Main Variables:
		- alphas increases to 0.9 from 0.8; it starts at the first update of a1n
	- Temperature:
		- Cold temperature decreasing - not as set in T1\_0 in main simulation.
		- temperature blows up towards the end of the simulation.
		- To ensure the positive part of the temperature, use abs( dB/dalpha ) instead of a signed value. Hack makes sense because plastic deformations work both ways.

- Try different boundary conditions (shock tube).

- Set rhos to constant artificially - answer: what does the simulation solve?
	- temperature gets too low
	- (ps - pg) is very large. Why?

- Fix solver of 4-equation model:
	- Use HLLC (same as original solver - check Toro)
		- doesn't work because because Fhllc needs to be changed (look for rho\*alphas)
	- How should temperature be computed without partial ps, es and rhos?
		- Find single-phase Hugoniot for original solver too.
	- Find Euler solver or write solver yourself. 


* Fundamental issues:
	- Fix overheating at boundary (bump due to different shock speeds?).
	- dt0 does not converge for large dx/dt
	- Revise artificially tweaked parameters: cv1, Ri0, Tc.

* Extend:
	- set n = 5 (optimization)
	- The temperature profile along pore is linear (extend to 5th order)

* Speed and Accuracy:
	- No need to open/close file if RAM is large (on cluster)
	- debug Newton code
	- Vectorize Fhllc
	- Extend to 2nd order (based on appendix).
	- make file reading faster.

* Read-write:
	- FIX savemode read and write, fa2a.m and fa2a2.m
	- Save and read reduced file (xx, alpha, p, T) skipping time steps

)
