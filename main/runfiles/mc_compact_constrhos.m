clear variables
close all
clc

% Simple Monte Carlo Simulation with initial porosity distribution.

pathfile
addpath(mainfunc_dir) % add path of helper functions
global gam1 gam2 pinf1 pinf2 c1 c2 rhos
initial_constrho


time 		= linspace(0, tend, tsteps);
xx          	= linspace(0, xend, ncells);
u0vec 		= zeros(1, length(xx)); u0vec(1) = u0;
p0vec		= p0 .* ones(1, length(xx));
varformat 	= '%8.5g';


% IF mc_count is set to 1, save as single_... in results/data/
if mc_count > 1
	newsim 		= 1;	% Whether to add a new folder for new MC set
	foldername0 = 'sim'; 
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
