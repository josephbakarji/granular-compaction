clear variables
close all
clc

% Simple Monte Carlo Simulation with initial porosity distribution.

pathfile
addpath(mainfunc_dir) % add path of helper functions
global gam1 gam2 pinf1 pinf2 c1 c2 dt0f_count rhos

dt0f_count = 0;

mc_count 	= 1;		% Number of MC realizations

tsteps 		= 1000; 	% number of steps
ncells 		= 100; 		% number of grid points
tend 		= 34e-6; 	% simulation time (s)
xend 		= 0.016; 	% domain size (m)

plotit 		= 1; 		% plot if 1.
plotmode 	= 0;		% plotmodes 0 or 1
plotfrq 	= 1;		% number of times steps per plot
printfrq 	= 100;		% number of times steps per print
savemode 	= 3;		% determines which variables to be saved: see save_result.m
savefile    	= 0;        % set to 1 to save data and to 0 to not save it

no_relax 	= 0;		% debug: set to 1 to test without relaxation step
compute_temp 	= 1;		% set to 1 to compute temperature distirbution
temp_mode	= 0;		% set to 0 for Saurel (2018), set to 1 for Carroll-Holt dissipation.
opt_tol 	= 1e-9;		% optimization (relaxation) tolerance
newton_opt 	= 0;		% use Newton's method (doesn't work)
bc          	= 1;         	% whether BCs are used (1 by default)
updatemesh 	= 1; 		% whether mesh is dynamically updated (default)

% heat dissipation variables
% see "Multiscale multiphase modeling of detonations in condensed energetic materials" (Saurel 2018)
muacc 		= 0.01;
Ls          	= 0.4538; 	% for pbx9501 %% use different name
n           	= 5;
as0 		= 0.5; % 0.98 in app. a (saurel 2018).
Ri0 		= 0.7e-5;

% Carroll-Holt dissipation model in spherical shell
% See "Dynamics of Heterogeneous Materials" (Nesterenko 2013) P. 153
etam            = 2e-3;
Y1          	= 400e6; 
B           	= 5765;
Tm          	= 1000; % ??? random

% problem params: for stiffened gas EOS
gam1 		= 5.5;
gam2 		= 1.4;
pinf1 		= 6e8; 	% pa % CHOSEN TO MAKE TEMPERATURE POSITIVE (HMX=31e8):
pinf2 		= 0;		% pa 
c1          	= 1500; 	% (m/s)  speed of sound in water !! aprox of hmx? !!
c2          	= 343;		% (m/s) speed of sound in air
cv1 		= 1444; % heat capacity = 1444.2 in App. A (Saurel 2018) - should change with temperature (assuming constant)
cv2		= 718; % J/K/kg (air - different for polymer binder)
T1_0		= 300; % Assuming room temperature
T2_0		= 300; % Assuming room temperature

m 		= 1.01;
np 		= 1.4;
a 		= 2000; % Tweaked

% Initial conditions
u0          	= 60;		% velocity in m/s
p0          	= 2e9;		% pressure in pa
rhos		= 1903
rho1_0 		= rhos;  	% hmx density in kg/m^3
rho2_0 		= 1;		% air density in kg/m^3

avg_alpha	= 0.8;		% porosity average
std_alpha	= 0.05;		% porosity std 
distribution	= 0;		% distribution type (see get_initial_por.m)

% ncells_mat = [60, 60, 120, 120, 120];
% tsteps_mat = [300, 600, 600, 1000, 1200];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for i = 1:length(ncells_mat)

%     ncells      = ncells_mat(i);
%     tsteps      = tsteps_mat(i);

time 		= linspace(0, tend, tsteps);
xx          	= linspace(0, xend, ncells);
u0vec 		= zeros(1, length(xx)); u0vec(1) = u0;
p0vec		= p0 .* ones(1, length(xx));
varformat = '%8.5g';

%% problem parameters

% IF mc_count is set to 1, save as single_... in results/data/
if mc_count > 1
	newsim 		= 1;	% Whether to add a new folder for new MC set
	foldername0 	= 'sim'; 
	savedir 	= [MonteCarlo_dir, foldername0];
	msgid 		= 'MATLAB:MKDIR:DirectoryExists';

	% Automatically add a directory for new set of simulations
	if newsim == 1
		k = 0;
		while(strcmp(msgid,'MATLAB:MKDIR:DirectoryExists'))
			k = k + 1;
			foldername = [foldername0, num2str(k)];
			savedir = [MonteCarlo_dir, foldername, '/'];
			[status, msg, msgid] = mkdir(savedir);
		end
	else
		[status, msg, msgid] = mkdir(savedir);
		if(status == 0)
			disp('folder could not be created')
		else
			disp(msg)
		end
	end
	nameprefix = 'MC';

else
	savedir = [data_dir, 'singles/'];
	nameprefix = 'single_';
end

disp('eooooooo')
for i = 1:mc_count

	if(i == 1)
		savemode = 1;
	else
		savemode = 3;
	end

	disp(['simulation n.: ', num2str(i)])
	alpha1_0 = get_initial_por(distribution, avg_alpha, std_alpha, xx);

	if savefile
		savefilename = makefilename(distribution, i, nameprefix);
	else
		savefilename = '';
	end

	Uo = make_uo(xx, alpha1_0, u0vec, p0vec, rho1_0, rho2_0, gam1, gam2, pinf1, pinf2, cv1, cv2, T1_0, T2_0);
	saur_main(savefilename, savedir, savemode, no_relax, compute_temp, temp_mode, plotit, plotfrq, printfrq, plotmode,...
		Uo, time, xx, bc, u0, updatemesh, varformat, opt_tol, newton_opt,...
		gam1, gam2, pinf1, pinf2, c1, c2, cv1, muacc, Ls, n, as0, Ri0, etam, Y1, B, Tm,...
		avg_alpha, std_alpha, distribution);

end
