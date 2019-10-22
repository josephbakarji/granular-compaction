clear all
close all
clc

%% Load Data
% Run get_2pntcorr.m if not working
load('time_corr_data.mat')

% Run mc_analysis_kde if not working
load('kde_analysis_Tstruct.mat')
load('kde_analysis_mfrac_str.mat')
load('kde_analysis_vars.mat')

%% 
addpath('../../runfiles/')
addpath('../../helpfunc/')
pathfile
foldernamelist = {'gauss_001_mf','gauss_002_mf', 'gauss_003_mf', 'gauss_004_mf', 'gauss_005_mf'};

sig_norm = std_a./avg_a;

%% Pre-processing

% Debug - NEGATIVE TEMPERATURES
[xi, yi, zi] = ind2sub(size(Tstruct{5}),find(Tstruct{5} <0));
for i = 1:length(Tstruct)
    [xi, yi, ri] = ind2sub(size(Tstruct{i}),find(Tstruct{i} <0));
    riu = unique(ri);
    idx = 1:size(Tstruct{i},3);
    Tnzidx{i} = setdiff(idx, riu);
end
%TT2 = Tstruct{5}(:, yi(56), zi(56));

% Variable Parameters
Tign = [600, 650, 700, 750];
sig_norm = std_a./avg_a; % Can be read from MC files

% CDF of temperature at a single point (xhs, ths)
xhs = 0.008599; xhsidx = find(xx_var{1}<=xhs, 1, 'last');
ths = 2.02e-5; thsidx = find(time_var{1}<=ths, 1, 'last');
disp('Computing CDF for point (xhs, ths)')
for j = 1:length(Tstruct)
    Temphs = reshape(Tstruct{j}(xhsidx, thsidx, Tnzidx{i}), 1, length(Tnzidx{i}));
    % Uses KDE
    [Tpdf_kde{j}, Tdis{j}] = ksdensity(Temphs, 'Function', 'pdf');
    [Tcdf_kde{j}, ~] = ksdensity(Temphs, Tdis{j}, 'Function', 'cdf');
    
end

% Probability of Ignition
P_ign = zeros(length(Tign), length(Tstruct));
for j = 1:length(Tstruct)
    for i = 1:length(Tign)
        TgTi = Tstruct{j}>Tign(i);
        P_ign(i,j) = nnz(TgTi)/numel(TgTi);
    end
end

% Average temperature
disp('Computing Average temperature (x, t)')
for i = 1:length(Tstruct)
    T_av{i} = mean(Tstruct{i}(:, :, Tnzidx{i}), 3);
end

% Average Rate of Reaction
Tdag0 = 2.65e4; %K
Z0 = 5e19; % (1/s)
Rrate = @(lams, Z, Tdag, Ts) (1 - lams) .* Z .* exp( - Tdag./Ts); % Reaction Rate Bdzil
for k = 1:length(Tstruct) % variances
	disp(['KDE map for $\sigma$ = ', num2str(sig_norm(k))])
	Rrate_tot = Rrate(mfrac_str{k}(:,:,Tnzidx{k}) , Z0, Tdag0, Tstruct{k}(:,:,Tnzidx{k}));
	tic;
	for i = 1:size(Tstruct{1}, 1) % x
		for j = 1:size(Tstruct{1}, 2) % time
            avg_rate{k}(i, j) = mean(Rrate_tot(i, j, :));
            avg_T{k}(i, j) = mean(Tstruct{k}(i, j, Tnzidx{k}));
		end
	end
	kdemap_time = toc;
	disp(['comptation time = ', num2str(kdemap_time)]);
end

% TOTAL average rate
for i = 1:length(avg_rate)
    totavg_rate(i) = mean(avg_rate{i}(:));
end

%% Plotting Options
xysize = 20;
legendsize = 18; 

%% Time correlation Figure


for j = 1:length(foldernamelist)
    [time_cdf{j}, ctx{j}] = ksdensity(Thist_nozero{j}, 'Function', 'cdf', 'Support', 'positive','BoundaryCorrection','reflection');
    [time_pdf{j}, tx{j}] = ksdensity(Thist_nozero{j}, 'Function', 'pdf', 'Support', 'positive','BoundaryCorrection','reflection');
end

fig0 = figure;
% for j = 1:length(foldernamelist)
%     plot(ctx{j}./1e-6, time_cdf{j}, 'LineWidth', 2)
%     hold on
%     sigleg{j} = ['$\bar \sigma_g = $', num2str(sig_norm(j))];
% end

plot(ctx{1}./1e-6, time_cdf{1}, '-', 'LineWidth', 3); hold on
plot(ctx{3}./1e-6, time_cdf{3}, '-.', 'LineWidth', 3); hold on
plot(ctx{5}./1e-6, time_cdf{5}, '--', 'LineWidth', 3);
ax = gca; 
ax.FontSize = 13;
sigleg{1} = ['$\bar \sigma_g = $ ', num2str(sig_norm(1))];
sigleg{2} = ['$\bar \sigma_g = $ ', num2str(sig_norm(3))];
sigleg{3} = ['$\bar \sigma_g = $ ', num2str(sig_norm(5))];

ylabel('$\mathcal F_t[\tau | T_\mathrm{ps}>T_\mathrm{ig}]$', 'Interpreter', 'Latex','FontSize', xysize);
xlabel('Time spent above threshold, $\tau$ ($\mu$s)', 'Interpreter', 'Latex','FontSize', xysize);
legend(sigleg, 'Interpreter', 'Latex', 'FontSize', legendsize, 'Location', 'SouthEast')
legend boxoff




fig0.PaperPosition = [0 0 6 9];
fig0.PaperSize = [6, 9];
orient(fig0,'landscape')
print(fig0,[plot_dir, '07_timecorr.pdf'],'-dpdf')

%% Initiation Probability with Variance

figv = figure;
subplot(1, 2, 1)
% for j = 1:length(Tstruct)
%     plot(Tdis{j}, Tcdf_kde{j}, 'LineWidth', 2)
%     hold on
%     legsig{j} = ['$\bar \sigma_g =$', num2str(sig_norm(j))];
% end

    plot(Tdis{1}, Tcdf_kde{1}, '-', 'LineWidth', 2, 'MarkerSize', 5); hold on
    plot(Tdis{3}, Tcdf_kde{3}, '-.', 'LineWidth', 2, 'MarkerSize', 4); hold on
    plot(Tdis{5}, Tcdf_kde{5}, '--', 'LineWidth', 2, 'MarkerSize', 3); hold on
    ax = gca; 
ax.FontSize = 13;

    legsig{1} = ['$\bar \sigma_g =$ ', num2str(sig_norm(1))];
    legsig{2} = ['$\bar \sigma_g =$ ', num2str(sig_norm(3))];
    legsig{3} = ['$\bar \sigma_g =$ ', num2str(sig_norm(5))];

xlabel('Sample Space Temperature, $\hat T$ (K)', 'Interpreter', 'Latex', 'FontSize', xysize)
ylabel('CDF, $\mathbf{P}[T_\mathrm{ps}(x_\mathrm{h}, t_\mathrm{h})<\hat T]$', 'Interpreter', 'Latex', 'FontSize', xysize)
legend(legsig, 'FontSize', legendsize, 'Interpreter', 'Latex', 'Location', 'Southeast')
legend boxoff

subplot(1, 2, 2)
% for i = 1:length(Tign)
%     plot(sig_norm, P_ign(i, :), '.-', 'MarkerSize', 30, 'LineWidth', 2); hold on;
%     
%     leg{i} = ['$T_\mathrm{ig} =$ ', num2str(Tign(i)), ' K'];
% end

plot(sig_norm, P_ign(1, :), '.-', 'MarkerSize', 30, 'LineWidth', 2); hold on;
plot(sig_norm, P_ign(2, :), '^-', 'MarkerSize', 8, 'LineWidth', 2); hold on;
plot(sig_norm, P_ign(3, :), 'o-', 'MarkerSize', 8, 'LineWidth', 2); hold on;
plot(sig_norm, P_ign(4, :), '*-', 'MarkerSize', 10, 'LineWidth', 2);

ax = gca; 
ax.FontSize = 13;

leg{1} = ['$T_\mathrm{ig} =$ ', num2str(Tign(1)), ' K'];
leg{2} = ['$T_\mathrm{ig} =$ ', num2str(Tign(2)), ' K'];
leg{3} = ['$T_\mathrm{ig} =$ ', num2str(Tign(3)), ' K'];
leg{4} = ['$T_\mathrm{ig} =$ ', num2str(Tign(4)), ' K'];

xlabel('Coefficient of Variation, $\bar \sigma_g$', 'Interpreter', 'Latex', 'FontSize', xysize)
ylabel('Ignition Probability, $\mathbf{P}[T_\mathrm{ps}(x_\mathrm{h}, t_\mathrm{h})>T_\mathrm{ig}]$', 'Interpreter', 'Latex', 'FontSize', xysize)
legend(leg, 'FontSize', legendsize, 'Interpreter', 'Latex', 'location', 'northwest')
legend boxoff

figv.PaperPosition = [0 0 6 14];
figv.PaperSize = [6, 14];
orient(figv,'landscape')
print(figv,[plot_dir, '06-var_cdf.pdf'],'-dpdf')


%% Total avg rate

figrt =figure;

subplot(1, 2, 1)
plot(sig_norm, totavg_rate, '.-', 'LineWidth', 3, 'MarkerSize', 30);
ax = gca; 
ax.FontSize = 13;


xlabel('Coefficient of Variation, $\overline \sigma_g$', 'Interpreter', 'Latex', 'FontSize', xysize)
ylabel('$\langle \dot R \rangle_s$', 'Interpreter', 'Latex', 'FontSize', xysize)


subplot(1, 2, 2)
% picks = [1, 3, 5];
% linestyles = ['-o', '.-', '^-'];
% i = 0;
% for j = picks
%     i = i + 1;
% 	plot(xx_var{j}.*100, T_av{j}(:, floor(length(time_var{j})/2)), '.-', 'LineWidth', 2, 'MarkerSize', 16)
% 	hold on
%     legsig{j} = ['$\bar \sigma_g = $ ', num2str(sig_norm(j))];
% end

plot(xx_var{1}.*100, T_av{1}(:, floor(length(time_var{1})/2)), '.-', 'LineWidth', 2, 'MarkerSize', 20); hold on
plot(xx_var{3}.*100, T_av{3}(:, floor(length(time_var{3})/2)), '^-', 'LineWidth', 2, 'MarkerSize', 5); hold on
plot(xx_var{5}.*100, T_av{5}(:, floor(length(time_var{5})/2)), '.--', 'LineWidth', 2, 'MarkerSize', 20);
ax = gca; 
ax.FontSize = 13;

legsig{1} = ['$\bar \sigma_g = $ ', num2str(sig_norm(1))];
legsig{2} = ['$\bar \sigma_g = $ ', num2str(sig_norm(3))];
legsig{3} = ['$\bar \sigma_g = $ ', num2str(sig_norm(5))];

legend(legsig, 'FontSize', legendsize, 'Interpreter', 'Latex', 'Location', 'Northeast')
x22 = xlabel('$x$ (cm)', 'Interpreter', 'Latex', 'FontSize', xysize);
y22 = ylabel('$\langle T_{\mathrm{ps}} \rangle_s$ (K)', 'Interpreter', 'Latex', 'FontSize', xysize);
legend boxoff


figrt.PaperPosition = [0 0 6 14];
figrt.PaperSize = [6, 14];
orient(figrt, 'landscape')
print(figrt, [plot_dir, '04-avg_TandR.pdf'],'-dpdf')

%% AVG Temperature and Reaction Rate

xind0 = 3;
j = 3;
[X, Y] = meshgrid(xx_var{j}(xind0:end), time_var{j});
gfig = figure;

subplot(1, 2, 1)
pcolor(X./1e-2, Y./1e-6, avg_rate{j}(xind0:end, :)');
ax = gca; 
ax.FontSize = 13;


shading interp
c = colorbar;
% title('$\langle \dot R \rangle$',  'Interpreter', 'Latex', 'FontSize', xysize)
xlabel('$x$ (cm)', 'Interpreter', 'Latex', 'FontSize', xysize)
ylabel('time ($\mu$s)', 'Interpreter', 'Latex', 'FontSize', xysize)
c.Label.String = '$\langle \dot R \rangle_s$';
c.Label.Interpreter = 'latex';
c.Label.FontSize = legendsize;

subplot(1, 2, 2)
pcolor(X./1e-2, Y./1e-6, avg_T{j}(xind0:end, :)');
ax = gca; 
ax.FontSize = 13;

shading interp
c2 = colorbar;
% title('$\langle T_\mathrm{ps} \rangle$',  'Interpreter', 'Latex', 'FontSize', xysize)
xlabel('$x$ (cm)', 'Interpreter', 'Latex', 'FontSize', xysize)
ylabel('time ($\mu$s)', 'Interpreter', 'Latex', 'FontSize', xysize)
c2.Label.String = '$\langle T_\mathrm{ps} \rangle_s$';
c2.Label.Interpreter = 'latex';
c2.Label.FontSize = legendsize;

gfig.PaperPosition = [0 0 6 14];
gfig.PaperSize = [6, 14];
orient(gfig,'landscape')
print(gfig,[plot_dir, '05-Temp-Rrate.pdf'],'-dpdf')

%% Realization
% Read

	pertend = 0.6; % percentage of simulation
	
	f1 = 'singles/SINGLE_2_simID_1.txt';
	f2 = 'singles/SINGLE022_2_simID_1.txt';

	[FA1, params1] = readtxtfile(f1, data_dir);
	[FA2, params2] = readtxtfile(f2, data_dir);

%% plot
	tsteps1 = params1(4);
	tsteps2 = params2(4);
	compute_temp1 = params2(14);
	compute_temp2 = params2(14);
	tt1 = floor(pertend * tsteps1);
	tt2 = floor(pertend * tsteps2);

	[A1, xx1, T1] = fa2a(FA1, tt1, compute_temp1);
	[A2, xx2, T2] = fa2a(FA2, tt2, compute_temp2);

	% assuming same for both
	gam1 = params1(8);
	gam2 = params1(9);
	pinf1 = params1(10);
	pinf2 = params1(11);

	p1_1 = pk(pinf1, gam1, rho1(A1), e1(A1));
	p2_1 = pk(pinf2, gam2, rho2(A1), e2(A1));
	p1_2 = pk(pinf1, gam1, rho1(A2), e1(A2));
	p2_2 = pk(pinf2, gam2, rho2(A2), e2(A2));

	idx0 = 70;
	a11 = a1(A1);
	a12 = a1(A2);
	u1 = u(A1);
	u2 = u(A2);
	xx1c = xx1(idx0:end) - xx1(idx0);
	xx2c = xx2(idx0:end) - xx2(idx0);

	fig2 = figure;
    c = lines(2);
	subplot(1, 2, 1);
	w1 = plot(xx1c.*100, a11(idx0:end), 'o', 'MarkerFaceColor', c(1, :)); hold on;
	w2 = plot(xx2c.*100, a12(idx0:end), '^', 'MarkerFaceColor', c(2, :));
	x1 = xlabel('$x$ (cm)');
	y1 = ylabel('Solid Volume Fraction, $\alpha_s$');
	l1 = legend({['$\bar \sigma_g = $ ' num2str(sig_norm(2))], ['$\bar \sigma_g = $ ' num2str(sig_norm(4))]}, 'Location', 'northeast');
    legend boxoff
    xlim([0, 1.35])
    
    ax = gca; 
ax.FontSize = 13;
        subplot(1, 2, 2);
	w3 = plot(xx1c.*100, T1(idx0:end), 'o', 'MarkerFaceColor', c(1, :)); hold on;
	w4 = plot(xx2c.*100, T2(idx0:end), '^', 'MarkerFaceColor', c(2 ,:));
	x2 = xlabel('$x$ (cm)');
	y2 = ylabel('Pore Surface Temperature, $T_\mathrm{ps}$ (K)');
	l2 = legend({['$\bar \sigma_g = $ ' num2str(sig_norm(2))], ['$\bar \sigma_g = $ ' num2str(sig_norm(4))]}, 'Location', 'northeast');
    legend boxoff
    xlim([0, 1.35])

    ax = gca; 
ax.FontSize = 13;

	set([w1, w2, w3, w4], 'MarkerSize', 4)
	set([x1, x2, y1, y2], 'Interpreter', 'Latex', 'Fontsize', xysize)
    set([l1, l2], 'Interpreter', 'Latex', 'Fontsize', legendsize)
%	fig2.PaperPosition = [0 0 5 13];

	fig2.PaperPosition = [0 0 6 14];
	fig2.PaperSize = [6, 14];
	orient(fig2,'landscape')
	print(fig2,[plot_dir, '03-alpha_temp.pdf'],'-dpdf')

