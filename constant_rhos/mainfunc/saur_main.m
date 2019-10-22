function [Un, Ti] = saur_main(savefilename, savedir, savemode, no_relax, compute_temp, temp_mode, plotit, plotfrq, printfrq, plotmode,...
		Uo, time, xx, bc, vbc, updatemesh, varformat, opt_tol, newton_opt,...
		gam1, gam2, pinf1, pinf2, c1, c2, cv1, muacc, Ls, n, as0, Ri0, etam, Y1, B, Tm,...
		avg_alpha, std_alpha, distribution)

	% Main solver function for 5-equation model described in "Simple and efficient ..." (Saurel 2009)
	% Inputs:
	% savefilename: file name where data of simulation is stored; string (e.g. 'file_example') see makefilename.m
	% savedir: file directory where data of simulation is stored; string (e.g. 'results/data/')
	% savemode: selects data storage option: see make_array.m for more details
	% no_relax: set to 1 if debugging relaxation step. Set to 0 for normal solver.
	% compute_temp: set to 1 for computing pore surface temperature. 0 otherwise.
	% temp_mode: set to 0 for Saurel (2018), set to 1 for Carroll-Holt dissipation.
	% plotit: set to 1 to plot simulation results on the go
	% plotfrq: number of time steps per plot.
	% printfrq: number of times steps per print (of simulation info).
	% plotmode: [1 or 2] plotting mode options: see plot_var.m
	% Uo: Initial array of initial conditions size=(7 * length(xx)) 
	% time: array of time
	% xx: array of gridpoints.
	% bc: set to 1 when boundary conditions are used
	% vbc: velocity boundary condition.
	% updatemesh: set to 1 if mesh is updated dynamically according to velocity distirbution.
	% varformat: string format of saved values (floating point and precision)
	% opt_tol: optimization (relaxation step) tolerance
	% newton_opt: set to 1 to use newton gradient descent (doesn't work)
	% gam1, gam2, pinf1, pinf2: pressure EOS variables (see reference)
	% c1, c2: wave speeds
	% cv1: heat capacity
	% muacc: accumulation viscosity (see Saurel 2018)
	% Ls: related to diffusivity (see Saurel 2018)
	% n: determined shape of diffusion in sphere (~5) (see Saurel 2018)
	% as0: porosity of granular medium without confined pressure.
	% Ri0: initial pore radius.
	% etam: melting viscovity in CH model (see Nesterenko 2013)
	% Y1: yield strength (see Nesterenko 2013) 
	% B: constant (see Nesterenko 2013) 
	% Tm: melting temperature (see Nesterenko 2013)
	% avg_alpha: average solid volume fraction
	% std_alpha: std of solid volume fraction
	% distirbution: initial porosity statistical distribution (see get_initial_por.m)

	% TODO:
	% FIX overheating at boundary.

	%	ch_temp = 0;		% Use Nesterenko 2013 for temperature calculation 
	%	saurel_temp = 1;	% Use Saurel 2018 for temperature calculation	

	pathfile
	addpath([maindir, 'helpfunc/']) % add path of helper functions

	% Variable allocation
	tsteps = length(time);
	ncells = length(xx);
	tend = time(end);	
	xend = xx(end);		
	dt = time(2) - time(1);                 % time increments
	u0 = Uo(3, 1)./(Uo(1, 1) + Uo(2, 1));   % Velocity BC
	T0 = 300.*ones(1, size(Uo,2)); 		% FIX!! 
	Ti = T0;

	% save data to file
	save_result(savefilename, savedir, savemode, Uo, xx, varformat, 0, compute_temp, temp_mode, T0, xend, ncells, tend, tsteps, no_relax, c1, c2, gam1, gam2, pinf1, pinf2, updatemesh,...
		bc, vbc, opt_tol, newton_opt, cv1, muacc, Ls, n, as0, Ri0, avg_alpha, std_alpha, distribution);


	%% Temperature computation
	if compute_temp == 1

		N2 = a2(Uo(:,1))./(4/3*pi.*Ri0.^3); % assumed constant!
		Ri = ( 3.*a2(Uo)./(N2*4*pi) ).^(1/3);
		Rim1 = Ri;
		vs0 = 1./rho1(Uo);
		ps0 = pk(pinf1, gam1, rho1(Uo), e1(Uo));
		es0 = e1(Uo);
	end


	%% MAIN LOOP
	Un = Uo;
	for t = 1:tsteps

		ap1 = a1(Uo) .* pk( pinf1, gam1, rho1(Uo), e1(Uo));
		ap2 = a2(Uo) .* pk( pinf2, gam2, rho2(Uo), e2(Uo));

		for j = 2:ncells-1
			%dx = xx(j) - xx(j-1); 		
			dx = (xx(j+1) - xx(j-1))/2;
			for i = 1:4
				Un(i,j) = Uo(i,j) - dt/dx * ( Fhllc(Uo(:,j), Uo(:,j+1), i) - Fhllc(Uo(:,j-1), Uo(:,j), i) ) ;
			end
			[Far, Smr] = Fhllc(Uo(:,j), Uo(:,j+1), 5); % flipped Smr and Sml
			[Fal, Sml] = Fhllc(Uo(:,j-1), Uo(:,j), 5);
			Fpe1r = Fhllc(Uo(:,j), Uo(:,j+1), 6);
			Fpe1l = Fhllc(Uo(:,j-1), Uo(:,j), 6);
			Fpe2r = Fhllc(Uo(:,j), Uo(:,j+1), 7);
			Fpe2l = Fhllc(Uo(:,j-1), Uo(:,j), 7);
			%Un(5,j) = Uo(5,j) - dt/dx * ( Far - Fal - Uo(5,j)*(Smr - Sml) );
			Un(5,j) = Un(1,j)./rho1(Uo(1,j)); % sets volume fraction as a function of initial density
			Un(6,j) = Uo(6,j) - dt/dx * ( Fpe1r - Fpe1l + ap1(j) * (Smr - Sml) );
			Un(7,j) = Uo(7,j) - dt/dx * ( Fpe2r - Fpe2l + ap2(j) * (Smr - Sml) );
		end

		%% Relaxation step

		p10 = pk(pinf1, gam1,  rho1(Un), e1(Un));  
		p20 = pk(pinf2, gam2,  rho2(Un), e2(Un));  
		p0 = p10.*a1(Un) + p20.*a2(Un);	
		v1_0 = 1./rho1(Un); 				 
		v2_0 = 1./rho2(Un); 			
		pI0 = pI(rho1(Un), c1, p10, rho2(Un), c2, p20);

		if no_relax == 1    % Debug: no relaxation <- no source terms.
			Un(8,:) = p0;
			pr = p0;
			a1n = a1(Un);
			a2n = a2(Un);
		else
			v1 = @(p, idx) vk_s(gam1, pinf1, p10(idx), v1_0(idx), p, pI0(idx));
			v2 = @(p, idx) vk_s(gam2, pinf2, p20(idx), v2_0(idx), p, pI0(idx));
			optf= @(p, idx) Un(1,idx) .* v1(p, idx) + Un(2,idx) .* v2(p, idx) - 1; 
			dv1 = @(p, idx) dvk_s(gam1, pinf1, p10(idx), v1_0(idx), p, pI0(idx));
			dv2 = @(p, idx) dvk_s(gam2, pinf2, p20(idx), v2_0(idx), p, pI0(idx));
			optdf= @(p, idx) Un(1,idx) .* dv1(p, idx) + Un(2,idx) .* dv2(p,idx); 
			options = optimset('TolX', opt_tol);


			% Optimization done for each grid point according to III.4 (Saurel)
			for i = 1:size(Un, 2) 
				if(newton_opt)
					[pr(i), ex] = newton_s(@(p) optf(p, i) , @(p) optdf(p, i), p0(i), opt_tol, 200);			
					if(ex > opt_tol)
						disp('did not converge - newton')
						%[pvec, exvec] = newton(@(p) optf(p, i) , @(p) optdf(p, i), p0(i), opt_tol, 200);
						%pp = linspace(p0(i) - 1e8, p0(i) + 1e8, 1000);
						%figure
						%plot(pp, optf(pp))
					end
				else
					try
						[pr(i), ~, exitflag] = fzero(@(p) optf(p, i) ,[0, 1e12], options); % make sure right bound always includes zero
					catch
						[pr(i), ~, exitflag] = fzero(@(p) optf(p, i) ,[0, 1e15], options);
						pp = linspace(0, p0(i) + 1e8, 10000);
						f2 = figure;
						plot(pp, optf(pp, i));
						keyboard
					end

					if exitflag ~= 1
						disp('did not converge')	
					end
				end
			end
			% Compute new porosities
			%a1n = Un(1,:) .* vk_s(gam1, pinf1, p10, v1_0, pr, pI0); % (rho1*a1) * v1
			a1n = Un(1,:) ./ rho1(Un); % (rho1*a1) * v1
			a2n = Un(2,:) .* vk_s(gam2, pinf2, p20, v2_0, pr, pI0); % (rho2*a2) * v2

			rho1_0 = rho1(Un(:,2));
			rho2_0 = rho2(Un(:,2));


			%% BC and mesh

			%% INITIALIZATION

			% Reinitialization

			if bc == 1
				a1n(1) = a1n(2);
				a2n(1) = a2n(2);
				Un(1,1) = rho1_0 * a1n(1);              % Neuman BC 
				Un(2,1) = rho2_0 * a2n(1);              % Neuman BC 
				Un(3,1) = (2*u0 - u(Un(:,2))) * rho(Un(:,1)); % Dirichlet Velocity BC (u*rho)
				rhoe = Un(4,:) - 0.5* Un(3,:).^2 ./ rho(Un); % rho*e = rho*E - 1/2 * (rho*u)^2/rho
				rhoe(1) = rhoe(2);    % Total internal energy conserved and neumann enforced.
				Un(4,1) = rhoe(1) + 0.5* rho(Un(:,1)) * u0.^2;

				pn = ( rhoe - (a1n.*gam1.*pinf1 ./ (gam1 - 1) + a2n.*gam2.*pinf2./(gam2 - 1)))./( a1n./(gam1 - 1) + a2n./(gam2 -1) ); % new pressure. Saurel II.4
			else
				rhoe = Un(4,:) - 0.5* Un(3,:).^2 ./ rho(Un); % rho*e = rho*E - 1/2 * (rho*u)^2/rho
				pn = ( rhoe - (a1n.*gam1.*pinf1 ./ (gam1 - 1) + a2n.*gam2.*pinf2./(gam2 - 1)))./( a1n./(gam1 - 1) + a2n./(gam2 -1) ); % new pressure. Saurel II.4
			end
			% Update internal energies.
			Un(6,:) = ( pn + pinf1*gam1.*a1n + pinf2*gam2.*a2n - (gam2-1).*rhoe)./((gam1-1)-(gam2-1));  % Use stiffened gas EOS, and Energy eqn
			Un(7,:) = rhoe - Un(6,:);
			Un(8,:) = pn;
			Un(5,:) = a1n;   

		end

		if bc == 1
			Un(3,1) = (2*u0 - u(Un(:,2))) * rho(Un(:,2)); % Dirichlet Velocity BC (u*rho)
			Un(4,1) = rhoe(1) + 0.5 * rho(Un(:,1)) * u0^2;
			Un(:,end) = Un(:,end-1);
		end

		% Update node positions
		if updatemesh == 1
			xx = xx + u(Un).*dt;
		end

		% Calculate Temperature:
		%Ti = temp_psurf(Un, p10, p20, es0, ps0, vs0, N2, cv1, muacc, Ls, as0, c1, c2);
		[Ti, Rim1] = temp_ps2(Un, Rim1, dt, es0, ps0, vs0, N2, cv1, muacc, Ls, as0, c1, c2);

		%keyboard
		% Save data and plot 
		save_result(savefilename, savedir, savemode, Un, xx, varformat, 1, compute_temp, temp_mode, Ti)
		plot_var(Un, plotit, plotfrq, plotmode, xx, t, tend, tsteps, Ti)
		print_time(t, tsteps, tend, printfrq)
		Uo = Un;
	end
end

