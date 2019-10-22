clear variables
close all
clc

% Read Monte Carlo simulation results, builds PDF and CDFs, plots results of initiation probability.

addpath('../../runfiles/')
pathfile
%foldernamelist = {'gauss_001','gauss_002', 'gauss_003', 'gauss_004', 'gauss_005'};
%foldernamelist = {'gauss_001_mf','gauss_002_mf', 'gauss_003_mf', 'gauss_004_mf', 'gauss_005_mf'};
foldernamelist = {'soren_sim1', 'soren_sim2', 'soren_sim3'};


read_file_flag = 1; % If 0, load .mat file

%% Read Temperature Data
numfold = length(foldernamelist);
avg_a = zeros(1, numfold);
std_a = zeros(1, numfold);
Tstruct = cell(1, numfold);
mfrac_str = cell(1, numfold);

if(read_file_flag)
    for j = 1:numfold
        full_file = [MonteCarlo_dir, foldernamelist{j}, '_full', '.mat'];
        full_folder = [MonteCarlo_dir, foldernamelist{j}, '/'];

        if(exist(full_folder, 'file'))
            if(exist(full_file, 'file'))
                load(full_file, 'Tarray', 'avg_alpha', 'std_alpha', 'xx', 'time', 'mfrac_s')       
            else
                [Tarray, avg_alpha, std_alpha, xx, time, mfrac_s] = mc_read_full(foldernamelist{j});        
            end
            Tstruct{j} = Tarray;
            mfrac_str{j} = mfrac_s;
            avg_a(j) = avg_alpha;
            std_a(j) = std_alpha;
            xx_var{j} = xx;
            time_var{j} = time;

            clear Tarray avg_alpha std_alpha xx0 time0 mcfrac_s;
        else
            disp(['Folder ', full_folder, '  Doesnt Exist'])
        end
    end
    save('kde_analysis_Tstruct.mat', 'Tstruct', '-v7.3')
    save('kde_analysis_mfrac_str.mat', 'mfrac_str', '-v7.3')
    save('kde_analysis_vars.mat', 'avg_a', 'std_a', 'xx_var', 'time_var')
else
    load('kde_analysis_Tstruct.mat')
    load('kde_analysis_mfrac_str.mat')
    load('kde_analysis_vars.mat')
end

disp('Done reading ...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Debug - NEGATIVE TEMPERATURES
% 
% [xi, yi, zi] = ind2sub(size(Tstruct{5}),find(Tstruct{5} <0));
% TT2 = Tstruct{5}(:, yi(56), zi(56));
% figure
% plot(reshape(TT2, [1, length(TT2)]))
% xlabel('x')
% ylabel('Temperature')


%% Variable Parameters

Tign = [600, 650, 700, 750];
sig_norm = std_a./avg_a; % Can be read from MC files
for i = 1:3
    legendsor{i} = [num2str(length(xx_var{i})), ' grid cells, ', num2str(length(time_var{i})) , ' time steps'];
end

%% CDF of temperature at a single point (xhs, ths)

xhs = 0.00455;
ths = 1e-5;

disp('Computing CDF for point (xhs, ths)')
for j = 1:length(Tstruct)
     xhsidx(j) = find(xx_var{j} <= xhs, 1, 'last');
     thsidx(j) = find(time_var{j} <= ths, 1, 'last');
    
    Temphs = reshape(Tstruct{j}(xhsidx(j), thsidx(j), :), 1, size(Tstruct{j}, 3));
    % Uses KDE
    [Tpdf_kde{j}, Tdis{j}] = ksdensity(Temphs, 'Function', 'pdf');
    [Tcdf_kde{j}, ~] = ksdensity(Temphs, Tdis{j}, 'Function', 'cdf');
    
end

% PLOTTING
fig1 = figure;
subplot(1, 2, 1)
for j = 1:length(Tstruct)
    plot(Tdis{j}, Tcdf_kde{j}, 'LineWidth', 2)
    hold on
    legsig{j} = ['$\bar \sigma^2_\xi =$', num2str(sig_norm(j))];
end
xlabel('Temperature, $T$', 'Interpreter', 'Latex', 'FontSize', 18)
ylabel('$\mathcal{P}(T_{ps}<T; x_h, t_h)$', 'Interpreter', 'Latex', 'FontSize', 18)
title('CDF of temperature at $(x_h, t_h)$', 'Interpreter', 'Latex', 'FontSize', 18)
legend(legendsor, 'FontSize', 18, 'Interpreter', 'Latex', 'Location', 'Southeast')
%axis([300, 770, 0, 1])

P_ign = zeros(length(Tign), length(Tstruct));
for j = 1:length(Tstruct)
    for i = 1:length(Tign)
        TgTi = Tstruct{j}>Tign(i);
        P_ign(i,j) = nnz(TgTi)/numel(TgTi);
    end
end


%% Average Temperature

disp('Computing Average temperature (x, t)')
for i = 1:length(Tstruct)
    T_av{i} = mean(Tstruct{i}, 3);
end


% PLOTTING

fig0 = figure;
for j = 1:length(Tstruct) 
	plot(xx_var{j}, T_av{j}(:, floor(length(time_var{j})/2)), '.-', 'LineWidth', 2, 'MarkerSize', 16)
	hold on
    legsig{j} = ['$\bar \sigma_\xi^2 =$ ', num2str(sig_norm(j))];
end
legend(legendsor, 'FontSize', 18, 'Interpreter', 'Latex', 'Location', 'Southeast')
t2 = title('Pore Surface Temperature ', 'Interpreter', 'Latex', 'FontSize', 18);
x2 = xlabel('x (in m)', 'Interpreter', 'Latex', 'FontSize', 18);
y2 = ylabel('$T_{ps}$ (in K)', 'Interpreter', 'Latex', 'FontSize', 18);


fig0.PaperPosition = [0 0 6 12];
fig0.PaperSize = [6, 9];
orient(fig0, 'landscape')
print(fig0, [plot_dir, 'mean_temp.pdf'],'-dpdf')



%% Average Rate of Reaction

Tdag0 = 2.65e4; %K
Z0 = 5e19; % (1/s)
lams0 = 0.8;  % wrong -> !!! solid mass fraction data has to be stored !!!
Rrate = @(lams, Z, Tdag, Ts) (1 - lams) .* Z .* exp( - Tdag./Ts); % Reaction Rate Bdzil


for k = 1:length(Tstruct) % variances
	disp(['KDE map for $\sigma$ = ', num2str(sig_norm(k))])
	Rrate_tot = Rrate(mfrac_str{k}, Z0, Tdag0, Tstruct{k});
	tic;
	for i = 1:size(Tstruct{k}, 1) % x
		for j = 1:size(Tstruct{k}, 2) % time
%  			Rrate_temp = reshape(Rrate_tot(i, j, :), 1, size(Rrate_tot, 3));
%            T_temp = reshape(T_tot(i, j, :), 1, size(T_tot, 3));
% 			[Rate_pdf_kde, Rate_dis] = ksdensity(Rrate_temp, 'Function', 'pdf', 'NumPoints', 100); % Necessary to use KDE?
%           avg_rate{k}(i, j) = trapz( Rate_dis, Rate_pdf_kde);
            
            avg_rate{k}(i, j) = mean(Rrate_tot(i, j, :));
            avg_T{k}(i, j) = mean(Tstruct{k}(i, j, :));
		end
	end
	kdemap_time = toc;
	disp(['comptation time = ', num2str(kdemap_time)]);
end


%% PLOTTING

xind0 = 10;
for j = 1:length(Tstruct)
    [X, Y] = meshgrid(xx_var{j}(xind0:end), time_var{j});
    g{j} = figure;
    
    subplot(1, 2, 1)
    pcolor(X./1e-3, Y./1e-6, avg_rate{j}(xind0:end, :)');
    shading interp
    colorbar
    title(['Average Reaction Rate for $\sigma_\xi^2 = $', num2str(sig_norm(j))], 'Interpreter', 'Latex', 'FontSize', 18)
    xlabel('x (cm)', 'Interpreter', 'Latex', 'FontSize', 18)
    ylabel('time ($\mu$s)', 'Interpreter', 'Latex', 'FontSize', 18)
    
    subplot(1, 2, 2)
    pcolor(X./1e-3, Y./1e-6, avg_T{j}(xind0:end, :)');
    shading interp
    colorbar
    title(['Average Temperature for $\sigma_\xi^2 = $', num2str(sig_norm(j))], 'Interpreter', 'Latex', 'FontSize', 18)
    xlabel('x (cm)', 'Interpreter', 'Latex', 'FontSize', 18)
    ylabel('time ($\mu$s)', 'Interpreter', 'Latex', 'FontSize', 18)
    
    g{j}.PaperPosition = [0 0 6 13];
    g{j}.PaperSize = [6, 13];
    orient(g{j},'landscape')
    print(g{j},[plot_dir, 'Temp-Rrate-',num2str(j),'.pdf'],'-dpdf')
end




%% Initiation Probability with Variance

figv = figure;
subplot(1, 2, 1)

for j = 1:length(Tstruct)
    plot(Tdis{j}, Tcdf_kde{j}, 'LineWidth', 2)
    hold on
    legsig{j} = ['$\bar \sigma^2_\xi =$', num2str(sig_norm(j))];
end
xlabel('Temperature, $T$', 'Interpreter', 'Latex', 'FontSize', 18)
ylabel('$\mathcal{P}(T_{ps}<T; x_h, t_h)$', 'Interpreter', 'Latex', 'FontSize', 18)
title('CDF of temperature at $(x_h, t_h)$', 'Interpreter', 'Latex', 'FontSize', 18)
legend(legendsor, 'FontSize', 14, 'Interpreter', 'Latex', 'Location', 'Southeast')


subplot(1, 2, 2)
for i = 1:length(Tign)
    plot(sig_norm, P_ign(i, :), '.-', 'MarkerSize', 30, 'LineWidth', 2);
    hold on;
    leg{i} = ['$T_{ig} =$ ', num2str(Tign(i)), ' K'];
end

title('Dependence on Variance', 'Interpreter', 'Latex', 'FontSize', 18)
xlabel('Variance, $\sigma_\xi^2$', 'Interpreter', 'Latex', 'FontSize', 18)
ylabel('Ignition Probability $\mathcal{P}(T_{ps}>T_{ig}; x_h, t_h)$', 'Interpreter', 'Latex', 'FontSize', 18)
legend(leg, 'FontSize', 14, 'Interpreter', 'Latex', 'location', 'northwest')


figv.PaperPosition = [0 0 6 13];
figv.PaperSize = [6, 13];
orient(figv,'landscape')
print(figv,[plot_dir, 'var_cdf_kde.pdf'],'-dpdf')



%% Tmperature CDF - PDF

fig2 = figure;
subplot(1, 2, 1)
for j = 1:length(Tstruct)
    plot(Tdis{j}, Tcdf_kde{j}, 'LineWidth', 2)
    hold on
end
xlabel('Temperature', 'Interpreter', 'Latex', 'FontSize', 18)
ylabel('CDF $\mathcal{P}(T<\tilde T; x_h, t_h)$', 'Interpreter', 'Latex', 'FontSize', 18)
title('CDF of temperature at $(x_h, t_h)$', 'Interpreter', 'Latex', 'FontSize', 18)
legend(legendsor, 'FontSize', 18, 'Interpreter', 'Latex', 'location', 'southeast')
%axis([300, 770, 0, 1])

subplot(1, 2, 2)
for j = 1:length(Tstruct)
	plot(Tdis{j}, Tpdf_kde{j}, 'LineWidth', 2)
    hold on
end
xlabel('Temperature', 'Interpreter', 'Latex', 'FontSize', 18)
ylabel('PDF $f(T; x_h, t_h)$', 'Interpreter', 'Latex', 'FontSize', 18)
title('PDF of temperature at $(x_h, t_h)$', 'Interpreter', 'Latex', 'FontSize', 18)
legend(legendsor, 'FontSize', 18, 'Interpreter', 'Latex')
%xlim([300, 750])



fig2.PaperPosition = [0 0 6 12];
fig2.PaperSize = [6, 12];
orient(fig2,'landscape')
print(fig2,[plot_dir, 'cdf_pdf_kde.pdf'],'-dpdf')

