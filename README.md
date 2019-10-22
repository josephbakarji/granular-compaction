
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


# Research Tasks #
## ACTIVE ##

To fix:
	* Main Variables:
		- alphas increases to 0.9 from 0.8; it starts at the first update of a1n
	* Temperature:
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


## Programming Tasks ##
* Publish on github

* Fundamental issues:
	- Fix overheating at boundary (bump due to different shock speeds?).
	- dt0 does not converge for large dx/dt
	- Revise artificially tweaked parameters: cv1, Ri0, Tc.

* Extend:
	- set n = 5 (optimization)
	- The temperature profile along pore is linear (extend to 5th order)
	- Implement Favrie 'configurational stress' terms (?)

* Speed and Accuracy:
	- No need to open/close file if RAM is large (on cluster)
	- debug Newton code
	- Vectorize Fhllc
	- Extend to 2nd order (based on appendix).
	- make file reading faster.

* Read-write:
	- FIX savemode read and write, fa2a.m and fa2a2.m
	- Save and read reduced file (xx, alpha, p, T) skipping time steps

---------
## Notes ##

- spt30:	
	- swtiched to temp\_ps2.m for temperature calculation
	- solid energy is now increasing near boundary.
	- make\_uo.m changed to compute partial pressures based on room temperature.
	- Take note: the stiffened gas EOS along with the es = Ts * Cvs heat capacity equation have the unknown (Ts, es, Cvs, ps, rhos, gams, pinfs). At the beginning we known gams, pinfs, rhos, Cvs(?). The question is whether we know ps or Ts at the beginning and which one we should set first. It makes sense to start with a room temperature Ts=300, then find the corresponding es and then ps; otherwise, if ps is taken as room pressure, the temperature will be very high for Cvs=1444 as given by Saurel (2018 apdx).
	- Parameters to be tweaked:
		- pinf1, cv1, p0, a (dBconf)
		- 6e8, 1444, 2e9, 5000 (works)

- oct3:
	- (old) Issue: Mismatch in multiscale model. Mesoscale heat model uses pressure disequilibrium while macroscale model is an equilibrium model.
		- solution: dr/dt was used instead of (ps - pg) to calculate the temperature
	
