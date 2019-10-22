clear variables
close all
clc

% Notes:
addpath('../../runfiles/')

pathfile
%foldernamelist = {'gauss_001','gauss_002', 'gauss_003', 'gauss_004', 'gauss_005'};
foldernamelist = {'gauss_001_mf','gauss_002_mf', 'gauss_003_mf', 'gauss_004_mf', 'gauss_005_mf'};

Tig = 650;

read_flag = 0;
if(read_flag)
    for j = 1:length(foldernamelist)
        full_folder = [MonteCarlo_dir, foldernamelist{j}, '/'];
        [Tsup_time, xx, time, avg_alpha, std_alpha] = mc_read_timecorr(foldernamelist{j}, Tig);
        Tstruct{j} = Tsup_time;
        Thist{j} = reshape(Tstruct{j}, [1, numel(Tstruct{j})]);
        Thist_nozero{j} = nonzeros(Thist{j});
        ratnonzero(j) = length(Thist_nozero{j})/length(Thist{j});
        avg_a(j) = avg_alpha;
        std_a(j) = std_alpha;
        disp(['ratio of grid points where temperature exceeds limit: ', num2str(ratnonzero(j))])
    end
    save('time_corr_data.mat')
else
    load('time_corr_data.mat')
end

sig_norm = std_a./avg_a;


%% Plot


for j = 1:length(foldernamelist)
    [time_cdf{j}, ctx{j}] = ksdensity(Thist_nozero{j}, 'Function', 'cdf', 'Support', 'positive','BoundaryCorrection','reflection');
    [time_pdf{j}, tx{j}] = ksdensity(Thist_nozero{j}, 'Function', 'pdf', 'Support', 'positive','BoundaryCorrection','reflection');
	%figure
	%title([foldernamelist{j}, 'for a simulation time: ', num2str(time(end))])
	%ylabel('number of occurences for all x');
	%xlabel('time spent above threshold');
end

fig0 = figure;
subplot(1, 2, 2)
for j = 1:length(foldernamelist)
    plot(tx{j}./1e-6, time_pdf{j}, 'LineWidth', 2)
    hold on
    sigleg{j} = ['$\bar \sigma^2_\xi =$', num2str(sig_norm(j))];
end
	title(['PDF - sim. time = ', num2str(time(end)./1e-6), ' $\mu$s'], 'Interpreter', 'Latex','FontSize', 18)
	ylabel('$P_t[\tau | T_{ps}>T_{ig}]$', 'Interpreter', 'Latex','FontSize', 18);
	xlabel('Time spent above threshold ($\mu$s)', 'Interpreter', 'Latex','FontSize', 18);
    legend(sigleg, 'Interpreter', 'Latex', 'FontSize', 16)
    
subplot(1, 2, 1)
for j = 1:length(foldernamelist)
    plot(ctx{j}./1e-6, time_cdf{j}, 'LineWidth', 2)
    hold on
    sigleg{j} = ['$\bar \sigma^2_\xi =$', num2str(sig_norm(j))];
end
	title(['CDF - sim. time = ', num2str(time(end)./1e-6), ' $\mu$s'], 'Interpreter', 'Latex','FontSize', 18)
	ylabel('$F_t[\tau | T_{ps}>T_{ig}]$', 'Interpreter', 'Latex','FontSize', 18);
	xlabel('Time spent above threshold ($\mu$s)', 'Interpreter', 'Latex','FontSize', 18);
    legend(sigleg, 'Interpreter', 'Latex', 'FontSize', 16, 'Location', 'SouthEast')
    
    
fig0.PaperPosition = [0 0 6 13];
fig0.PaperSize = [6, 13];
orient(fig0,'landscape')
print(fig0,[plot_dir, 'time_pdf.pdf'],'-dpdf')

% 
% numbins = 30;
% 
% for j = 1:length(foldernamelist)
% 	figure
% 	hh{j} = histogram(Thist_nozero{j}, numbins);
% 	title([foldernamelist{j}, 'for a simulation time: ', num2str(time(end))])
% 	ylabel('number of occurences for all x');
% 	xlabel('time spent above threshold');
% end
% 
% figure
% for j = 1:length(foldernamelist)
%     xv = hh{j}.BinEdges(1:end-1) + hh{j}.BinWidth/2;
%     yv = hh{j}.Values;
%     plot(xv, yv, 'LineWidth', 2)
%     hold on
% end
% 	title(['For a simulation time: ', num2str(time(end))])
% 	ylabel('number of occurences for all x');
% 	xlabel('time spent above threshold');
%     legend('std = 0.01', 'std = 0.02', 'std = 0.03', 'std = 0.04')
