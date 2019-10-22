clear variables
close all
clc

% Simple Monte Carlo Simulation with initial porosity distribution.

addpath(genpath('../mainfunc/')) % add path of helper functions
global gam1 gam2 pinf1 pinf2 c1 c2

mc_count 	= 1;		% Number of MC realizations

dt          = 1.25e-8;
%tsteps 		= 2000; 	% number of steps
ncells 		= 200; 		% number of grid points
tend 		= 1e-5; 	% simulation time (s) - % for 50 m/s
xend 		= 0.016; 	% domain size (m)
tsteps      = floor(tend ./ dt);


plotit 		= 1; 		% plot if 1.
plotmode 	= 0;		% plotmodes 0 or 1
plotfrq 	= 200;		% number of times steps per plot
printfrq 	= 500;		% number of times steps per print
savemode 	= 2;		% determines which variables to be saved: see save_result.m

no_relax 	= 0;		% debug: set to 1 to test without relaxation step
compute_temp= 1;		% set to 1 to compute temperature distirbution
opt_tol 	= 1e-7;		% optimization (relaxation) tolerance
newton_opt 	= 0;		% use Newton's method (doesn't work)
bc          = 1;         	% whether BCs are used (1 by default)
updatemesh 	= 1; 		% whether mesh is dynamically updated (default)

% heat dissipation variables
% see "Multiscale multiphase modeling of detonations in condensed energetic materials" (Saurel 2018)
cv1 		= 6600; 	% heat capacity = 1444.2 in App. A (Saurel 2018). 
muacc 		= 0.01;
as0 		= 0.5; % 0.98 in app. a (saurel 2018).
Ri0 		= 0.7e-5;
n 		= 5;



% Carroll-Holt dissipation model in spherical shell
% See "Dynamics of Heterogeneous Materials" (Nesterenko 2013) P. 153
etam		= 2e-3;
Y1 		= 400e6; 
B 		= 5765;
Tm		= 1000; % ??? random

% problem params: for stiffened gas EOS
gam1 	= 5.5;
gam2 	= 1.4;
pinf1 	= 31e8; 	% pa
pinf2 	= 0;		% pa 
c1 		= 1500; 	% (m/s)  speed of sound in water !! aprox of hmx? !!
c2 		= 343;		% (m/s) speed of sound in air

% Initial conditions
p0          = 1e6;		% pressure in pa
rho1_0 		= 1903;  	% hmx density in kg/m^3
rho2_0 		= 1;		% air density in kg/m^3



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time 		= linspace(0, tend, tsteps);
xx          = linspace(0, xend, ncells);
p0vec 		= p0 .* ones(1, length(xx));


varformat = '%8.5g';

%% problem parameters
%
	datadir 	= '../results/data/';
	foldername0 	= 'parameterexp1/';
	savedir 	= [datadir, foldername0];

	% Automatically add a directory for new set of simulations
	[status, msg, msgid] = mkdir(savedir);

	if(status == 0)
		disp('folder could not be created')
	else
		disp(msg)
	end

Ls_mat 		= linspace(0.35, 0.55, 6); %0.4538; 	% for pbx9501 %% use different name
u0_mat 		= [1, 10, 40, 70, 100];		% velocity in m/s
tend_vec    = [1e-4, 1e-4, 4e-5, 2e-5, 1.8e-5]; 
avg_alpha_mat	= [0.5, 0.7, 0.8, 0.9, 0.98];		% porosity average
std_alpha_mat	= [0.01, 0.02, 0.03, 0.05, 0.06];		% porosity std  
distribution= 1;		% distribution type (see get_initial_por.m)


for i = 1:length(Ls_mat)
    for j = 1:length(u0_mat)     
        
        
        
		u0 = u0_mat(j);
		Ls = Ls_mat(i);
        
        tend        = tend_vec(j);
        tsteps      = floor(tend/dt);
        time 		= linspace(0, tend, tsteps);
        
		avg_alpha = avg_alpha_mat(3);
		std_alpha = std_alpha_mat(2);

        disp(['u0 = ', num2str(u0), '; Ds = ', num2str(Ls)])
        
		u0vec	= zeros(1, length(xx)); 
		u0vec(1) = u0;

		disp(['simulation n.: ', num2str(i)])
		alpha1_0 = get_initial_por(distribution, avg_alpha, std_alpha, xx);

		savefilename = makefilename(distribution, (i-1)*length(u0_mat) + j, 'pe');


		Uo = make_uo(alpha1_0, p0vec, u0vec, rho1_0, rho2_0, gam1, gam2, pinf1, pinf2);
		[UN{i, j}, Ti{i, j}] = saur_main(savefilename, savedir, savemode, no_relax, compute_temp, plotit, plotfrq, printfrq, plotmode,...
			Uo, time, xx, bc, updatemesh, varformat, opt_tol, newton_opt,...
			gam1, gam2, pinf1, pinf2, c1, c2, cv1, muacc, Ls, n, as0, Ri0, etam, Y1, B, Tm,...
			avg_alpha, std_alpha, distribution);
        
        %keyboard

	end	
end

keyboard
figure
for j = 1:length(u0_mat)
    figure
	for i = 1:length(Ls_mat)
		plot_var(UN{i,j}, plotit, plotfrq, plotmode, xx, nan, tend, tsteps, Ti{i, j})
		hold on
	end
end


