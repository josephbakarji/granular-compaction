clear variables
close all
clc

% Single simulation with uniform initial porosity.

pathfile
addpath(mainfunc_dir) % add path of helper functions
global gam1 gam2 pinf1 pinf2 c1 c2 dt0f_count
initial_simple

time 		= linspace(0, tend, tsteps);
xx          	= linspace(0, xend, ncells);
u0vec 		= zeros(1, length(xx)); u0vec(1) = u0;
p0vec		= p0 .* ones(1, length(xx));
varformat 	= '%8.5g'; % format of saved data

savedir 	= [data_dir, 'singles/'];
nameprefix 	= 'single_';


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

Uo = make_uo(xx, alpha1_0, u0vec, p0vec, rho1_0, rho2_0, gam1, gam2, pinf1, pinf2);
solve1D(savefilename, savedir, savemode, compute_temp, temp_mode, plotit, plotfrq, printfrq, plotmode,...
	Uo, time, xx, bc, vbc, updatemesh, varformat, opt_tol, gam1, gam2, pinf1, pinf2, c1, c2, cv1,...
	muacc, Ls, n, as0, Ri0, avg_alpha, std_alpha, distribution)
