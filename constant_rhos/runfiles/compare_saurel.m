clear variables
close all
clc

% Validation test cases in "Simple and efficient relaxation methods for interfaces separating compressible fluids..." (Saurel 2009) 
% Cases:
% 1 - shock tube: unequal pressures, equal initial volume fraction
% 2 - shock tube: unequal pressures, unequal initial volume fraction


addpath(genpath('../mainfunc/')) % add path of helper functions

global gam1 gam2 pinf1 pinf2 c1 c2


plotit 		= 1; 	% plot if 1.
plotmode	= 1;
plotfrq 	= 1;
printfrq 	= 20;
savemode 	= 2;

tend 		= 250e-6;
xend 		= 1;
no_relax 	= 0;	% DEBUG: test without relaxation step if 1
compute_temp 	= 0; 
opt_tol 	= 1e-7;
newton_opt 	= 0;
bc 		= 0;    % Whether BCs are used (default)
updatemesh 	= 0; 	% Whether mesh is dynamically updated (default)

cv1 		= 4000; % = 1444.2 from app. A (Saurel 2018). 
muacc 		= 0.01;    
Ls 		= 0.4538; % For PBX9501 %% USE DIFFERENT NAME
n 		= 5;    
as0 		= 0.5; % 0.98 in app. A (Saurel 2018).
Ri0 		= 0.5e-5;

rho1_0 		= 1000;  % water density in kg/m^3
rho2_0 		= 1;    % air density in kg/m^3
gam1 		= 4.4;
gam2 		= 1.4;
pinf1 		= 6e8; % Pa
pinf2 		= 0; % kPa 
c1 		= 1500; % (m/s)  Speed of sound in water 
c2 		= 343; % (m/s) speed of sound in air

etam		= 2e-3;
Y1 		= 400e6; 
B 		= 5765;
Tm		= 1000; % ??? random
avg_alpha	= 0.8;       % porosity 
std_alpha	= 0.02;
distribution	= 0;

varformat = '%9.7g';
savedir = 'results/data/';

%% Problem Parameters

timexx = {[700, 7000], [300, 3000]};
casepre = {'finegrid', 'coarsegrid'};
casenum = 2;
dis_opts = 2;

% Simulations

k = 0;
for j = casenum 
	for i = dis_opts 
		k = k + 1; 
		disp(['case: ', num2str(j), ' , discretization: ', num2str(i)])

		savefilename = '';%	makefilename(j, k, ['compare_saur_pdiffn_', casepre{i}]);
		tsteps = timexx{i}(2);
		ncells = timexx{i}(1);
		time = linspace(0, tend, tsteps);
		xx = linspace(0, xend, ncells);

		if(j == 1) 
		%% 2 - shock-tube: unequal pressures, equal initial volume fraction
			xs = 0.5;
			u0= 0;  % m/s
			p_left = 1e8;
			p_right = 1e5;
			alpha1_0 = 0.5 .* ones(1, length(xx));

			u0vec = u0.*ones(1, length(xx)); 
			p0vec = zeros(1, length(xx));

			[~, sepind] = min(abs(xx - xs));
			p0vec(1:sepind) = p_left;	% Left pressure
			p0vec(sepind+1:end) = p_right;	% Right pressure


		elseif(j == 2)
		%% 3 - shock tube: unequal pressures, unequal initial volume fraction
			xs = 0.5;
			u0= 0;  % m/s
			p_left = 1e8;
			p_right = 1e5;
			alpha1_0 = 0.5 .* ones(1, length(xx));
			eps = 0.01;

			u0vec = u0.*ones(1, length(xx)); 
			p0vec = zeros(1, length(xx));
			[~, sepind] = min(abs(xx - xs));
			p0vec(1:sepind) = p_left;		% Left pressure
			p0vec(sepind+1:end) = p_right;	% Right pressure
			alpha1_0(1:sepind) = eps;
			alpha1_0(sepind + 1: end) = 1 - eps;
		end

		Uo = make_uo(alpha1_0, p0vec, u0vec, rho1_0, rho2_0, gam1, gam2, pinf1, pinf2);
		saur_main(savefilename, savedir, savemode, no_relax, compute_temp, plotit, plotfrq, printfrq, plotmode,...
			Uo, time, xx, bc, updatemesh, varformat, opt_tol, newton_opt,...
			gam1, gam2, pinf1, pinf2, c1, c2, cv1, muacc, Ls, n, as0, Ri0, etam, Y1, B, Tm,...
			avg_alpha, std_alpha, distribution)
	end
end
