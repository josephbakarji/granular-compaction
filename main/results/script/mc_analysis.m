clear variables
close all
clc

% Read Monte Carlo simulation results, builds PDF and CDFs, plots results of initiation probability.

addpath('../../runfiles/')

pathfile
foldernamelist = {'gauss_001','gauss_002', 'gauss_003', 'gauss_004'};


for j = 1:length(foldernamelist)
    full_file = [MonteCarlo_dir, foldernamelist{j}, '.mat'];
    full_folder = [MonteCarlo_dir, foldernamelist{j}, '/'];
    
    if(exist(full_folder, 'file'))
        if(exist(full_file, 'file'))
            load(full_file, 'histmat', 'Tdis', 'xx', 'time')       
        else
            [histmat, Tdis, xx, time] = mc_read(foldernamelist{j});        
        end
        Hstruct{j} = histmat;
        clear histmat;
    else
        disp(['Folder ', full_folder, '  Doesnt Exist'])
    end
end

%%%%%%%%%% READ single %%%%%%%%%%%%%%%

%[FA, params] = readtxtfile('single0_mean_0_simID_1.txt', '../data/');
%[A, xx, Tmean] = fa2a2(FA, floor(length(time)/2), 2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% PLOTS
Tign = 700;
Tign2 = 650;

%% 1 - CDF of temperature at a given time ths and position xhs
xhs = 0.008599; xhsidx = find(xx<=xhs, 1, 'last');
ths = 2.02e-5; thsidx = find(time<=ths, 1, 'last');

for j = 1:length(Hstruct)
    Temphs{j} = reshape(Hstruct{j}(xhsidx, thsidx, :), 1, length(Tdis));
    Temphstot{j} = sum(Temphs{j});
    Tcdf{j} = zeros(1, length(Temphs{j}));
    Tcdf{j}(1) = Temphs{j}(1);
    for i = 2:length(Temphs{j})
        Tcdf{j}(i) = Temphs{j}(i)/Temphstot{j} + Tcdf{j}(i-1);
    end
end

% Calculate PDF
cutshort = 9;
for i = 1:length(Hstruct)
	Tcdf_s{i} = smooth(Tcdf{i}, 17);
	for j = 1:length(Tcdf{i})-cutshort
		Tpdf{i}(j) = (Tcdf_s{i}(j+1) - Tcdf_s{i}(j))/(Tdis(j+1) - Tdis(j));
	end
	Tpdf_s{i} = smooth(Tpdf{i}, 2);
end
Tdis_df = Tdis(1:end-cutshort);


% Calculate Average Temperature

for i = 1:length(Hstruct)
	for j = 1:size(Hstruct{i}, 1)
		for k = 1:size(Hstruct{i}, 2)
			sumTtot = sum(Hstruct{i}(j,k,:));
			Ttotal = 0;
			for l = 1:size(Hstruct{i}, 3)
				Ttotal = Hstruct{i}(j,k,l) * Tdis(l) + Ttotal;
			end
			T_av{i}(j,k) = Ttotal./sumTtot;
		end
	end
end



% Plot Mean temperature as a function of x, given time, different init. porosity std.
fig0 = figure;
for j = 1:length(Hstruct) 
	plot(xx, T_av{j}(:, floor(length(time)/2)), '.-', 'LineWidth', 2, 'MarkerSize', 16)
	hold on
end
legend({'$\sigma_\xi^2 = 0.01$','$\sigma_\xi^2 = 0.02$','$\sigma_\xi^2 = 0.03$','$\sigma_\xi^2 = 0.04$'}, 'FontSize', 18, 'Interpreter', 'Latex', 'Location', 'Southeast')
t2 = title('Pore Surface Temperature ', 'Interpreter', 'Latex', 'FontSize', 18);
x2 = xlabel('x (in m)', 'Interpreter', 'Latex', 'FontSize', 18);
y2 = ylabel('$T_{ps}$ (in K)', 'Interpreter', 'Latex', 'FontSize', 18);




fig1 = figure;
subplot(1, 2, 1)
for j = 1:length(Hstruct)
    plot(Tdis, Tcdf_s{j}, 'LineWidth', 2)
    hold on
end
xlabel('Temperature, $T$', 'Interpreter', 'Latex', 'FontSize', 18)
ylabel('$\mathcal{P}(T_{ps}<T; x_h, t_h)$', 'Interpreter', 'Latex', 'FontSize', 18)
title('CDF of temperature at $(x_h, t_h)$', 'Interpreter', 'Latex', 'FontSize', 18)
legend({'$\sigma_\xi^2 = 0.01$','$\sigma_\xi^2 = 0.02$','$\sigma_\xi^2 = 0.03$','$\sigma_\xi^2 = 0.04$'}, 'FontSize', 18, 'Interpreter', 'Latex', 'Location', 'Southeast')
%axis([300, 770, 0, 1])


for j = 1:length(Hstruct)
    Tignidx = find(Tdis< Tign, 1, 'last');
    TotTign{j} = sum(Hstruct{j}(:, :, Tignidx:end), 3);
    TotT{j} = sum(Hstruct{j}, 3);
    IgnPr{j} = TotTign{j} ./ TotT{j};

    Tignidx2 = find(Tdis< Tign2, 1, 'last');
    TotTign2{j} = sum(Hstruct{j}(:, :, Tignidx2:end), 3);
    IgnPr2{j} = TotTign2{j} ./ TotT{j};
end


for i = 1:length(Hstruct)
    pign(i) = 1 - Tcdf{i}(Tignidx);
    pign2(i) = 1 - Tcdf{i}(Tignidx2);
end

variances = [0.01, 0.02, 0.03, 0.04];

subplot(1, 2, 2)
plot(variances, pign, '.-', 'MarkerSize', 30, 'LineWidth', 2);hold on;
plot(variances, pign2, '.-', 'MarkerSize', 30, 'LineWidth', 2)
title('Dependence on Variance', 'Interpreter', 'Latex', 'FontSize', 18)
xlabel('Variance, $\sigma_\xi^2$', 'Interpreter', 'Latex', 'FontSize', 18)
ylabel('Ignition Probability $\mathcal{P}(T_{ps}>T_{ig}; x_h, t_h)$', 'Interpreter', 'Latex', 'FontSize', 18)
legend({['$T_{ig} =$ ', num2str(Tign), ' K'], ['$T_{ig} =$ ', num2str(Tign2), ' K']}, 'FontSize', 16, 'Interpreter', 'Latex', 'location', 'northwest')


fig1.PaperPosition = [0 0 6 12];
fig1.PaperSize = [6, 12];
orient(fig1,'landscape')
print(fig1,[plot_dir, 'var_depend.pdf'],'-dpdf')

%%%
%%%
%%%

fig2 = figure;
subplot(1, 2, 1)
for j = 1:length(Hstruct)
    plot(Tdis, Tcdf_s{j}, 'LineWidth', 2)
    hold on
end
xlabel('Temperature', 'Interpreter', 'Latex', 'FontSize', 18)
ylabel('CDF $\mathcal{P}(T<\tilde T; x_h, t_h)$', 'Interpreter', 'Latex', 'FontSize', 18)
title('CDF of temperature at $(x_h, t_h)$', 'Interpreter', 'Latex', 'FontSize', 18)
legend({'$\sigma_\xi^2 = 0.01$','$\sigma_\xi^2 = 0.02$','$\sigma_\xi^2 = 0.03$','$\sigma_\xi^2 = 0.04$'}, 'FontSize', 18, 'Interpreter', 'Latex', 'location', 'southeast')
%axis([300, 770, 0, 1])

subplot(1, 2, 2)
for j = 1:length(Hstruct)
	plot(Tdis_df, Tpdf_s{j}, 'LineWidth', 2)
    hold on
end
xlabel('Temperature', 'Interpreter', 'Latex', 'FontSize', 18)
ylabel('PDF $f(T; x_h, t_h)$', 'Interpreter', 'Latex', 'FontSize', 18)
title('PDF of temperature at $(x_h, t_h)$', 'Interpreter', 'Latex', 'FontSize', 18)
legend({'$\sigma_\xi^2 = 0.01$','$\sigma_\xi^2 = 0.02$','$\sigma_\xi^2 = 0.03$','$\sigma_\xi^2 = 0.04$'}, 'FontSize', 18, 'Interpreter', 'Latex')
%xlim([300, 750])



fig2.PaperPosition = [0 0 6 12];
fig2.PaperSize = [6, 12];
orient(fig2,'landscape')
print(fig2,[plot_dir, 'cdf_pdf.pdf'],'-dpdf')

%% 2- CDF of ignition at a specific position x in time GIVEN that it initiates


%xhs = 0.008599; xhsidx = find(xx<=xhs, 1, 'last');
%
%for j = 1:length(Hstruct)
%    Temphs_t{j} = reshape( sum(Hstruct{j}(xhsidx, :, Tignidx:end), 3)./sum(Hstruct{j}(xhsidx, :, :), 3) , 1, length(time));
%    Temphstot_t{j} = sum(Temphs_t{j});
%    Tcdf_t{j} = zeros(1, length(time));
%    Tcdf_t{j}(1) = Temphs_t{j}(1);
%    for i = 2:length(Temphs_t{j})
%        Tcdf_t{j}(i) = Temphs_t{j}(i)/Temphstot_t{j} + Tcdf_t{j}(i-1);
%    end
%end 

%subplot(2, 2, 2)
%
%for j = 1:length(Hstruct)
%    plot(time, Tcdf_t{j}, 'LineWidth', 2)
%    hold on;
%end
%
%xlabel('time', 'Interpreter', 'Latex')
%ylabel('CDF $P(t<\bar{t} | T>T_{ig})$', 'Interpreter', 'Latex')
%title('CDF of ignition at a given position x', 'Interpreter', 'Latex')
%legend('\sigma_xi^2 = 0.01','\sigma_xi^2 = 0.02','\sigma_xi^2 = 0.03','\sigma_xi^2 = 0.04')
%
%
%
%%% 3 - CDF in x at a given time.
%
%% CDF of ignition at a specific time in x, GIVEN that it initiates
%ths = 2.02e-5; 
%thsidx = find(time<=ths, 1, 'last');
%
%for j = 1:length(Hstruct)
%    Temphs_x{j} = reshape( sum(Hstruct{j}(:, thsidx, Tignidx:end), 3)./sum(Hstruct{j}(:, thsidx, :), 3) , 1, length(xx));
%    Temphstot_x{j} = sum(Temphs_x{j});
%    Tcdf_x{j} = zeros(1, length(xx));
%    Tcdf_x{j}(1) = Temphs_x{j}(1);
%    for i = 2:length(Temphs_x{j})
%        Tcdf_x{j}(i) = Temphs_x{j}(i)/Temphstot_x{j} + Tcdf_x{j}(i-1);
%    end
%end 
% 
%subplot(2, 2, 3)
%
%for j = 1:length(Hstruct)
%    plot(xx, Tcdf_x{j}, 'LineWidth', 2)
%    hold on;
%end
%
%xlabel('time', 'Interpreter', 'Latex')
%ylabel('CDF $P(x<\bar{x} | T>T_{ig})$', 'Interpreter', 'Latex')
%title('CDF of ignition at a given time t', 'Interpreter', 'Latex')
%legend('\sigma_xi^2 = 0.01','\sigma_xi^2 = 0.02','\sigma_xi^2 = 0.03','\sigma_xi^2 = 0.04')
%

%% Probability of ignition as a function of variance


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Surface Probability of initiation

xidx = [0.2, 0.4];
tidx = [0.2, 0.5];
xidxrange = floor(xidx(1) * length(xx)) : floor(xidx(2) * length(xx));
tidxrange = floor(tidx(1) * length(time)) : floor(tidx(2) * length(time));
xnew = xx(xidxrange);
tnew = time(tidxrange);
IgnPrnew = IgnPr{3}(xidxrange, tidxrange);


% Cumulative
IgnCum = zeros(size(IgnPrnew));
IgnCum(1, 1) = IgnPrnew(1, 1);
for i = 1:length(xnew)
	for j = 1:length(tnew)
		if(i == 1)
			xcum = 0;
		else
			xcum = IgnCum(i - 1, j);
		end

		if(j == 1)
			tcum = 0;
		else
			tcum = IgnCum(i , j - 1);
		end

		IgnCum(i, j) = IgnPrnew(i, j) + xcum + tcum;
	end
end
IgnCum = IgnCum./IgnCum(end, end);



xnewint = linspace(xnew(1), xnew(end), 500);
tnewint = linspace(tnew(1), tnew(end), 500);
[X, Y] = meshgrid(xnew, tnew);
[Xq, Yq] = meshgrid(xnewint, tnewint);
Vq = interp2(X, Y, IgnPrnew', Xq, Yq);
Vq_smooth = imfilter(Vq,fspecial('average',[30 30])); 

fig3 = figure;
subplot(1, 2, 1)
surf(Xq, Yq, Vq_smooth)
%colorbar
shading interp
title('Probability of Initiation', 'Interpreter', 'Latex', 'FontSize', 16)
xlabel('x (in m)', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel('time (in s)', 'Interpreter', 'Latex', 'FontSize', 14)
%, $T_{ig}=$ ', num2str(Tign), ' K']


subplot(1, 2, 2)
pcolor(Xq, Yq, Vq)
colorbar
shading interp

title('Probability of Initiation', 'Interpreter', 'Latex', 'FontSize', 16)
xlabel('x (in m)', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel('time (in s)', 'Interpreter', 'Latex', 'FontSize', 14)


fig3.PaperPosition = [0 0 6 13];
fig3.PaperSize = [6, 13];
orient(fig3,'landscape')
print(fig3,[plot_dir, 'surface_prob.pdf'],'-dpdf')



fig4 = figure;

subplot(1, 2, 1)
for j = 1:length(Hstruct)
    plot(Tdis, Tcdf_s{j}, 'LineWidth', 2)
    hold on
end
xlabel('Temperature', 'Interpreter', 'Latex', 'FontSize', 18)
ylabel('CDF $\mathcal{P}(T<\tilde T; x_h, t_h)$', 'Interpreter', 'Latex', 'FontSize', 18)
title('CDF of temperature at $(x_h, t_h)$', 'Interpreter', 'Latex', 'FontSize', 18)
legend({'$\sigma_\xi^2 = 0.01$','$\sigma_\xi^2 = 0.02$','$\sigma_\xi^2 = 0.03$','$\sigma_\xi^2 = 0.04$'}, 'FontSize', 18, 'Interpreter', 'Latex', 'Location', 'Southeast')
%axis([300, 770, 0, 1])


subplot(1, 2, 2)
pcolor(Xq, Yq, Vq)
colorbar
shading interp

title(['Ignition Probability with $T_{ig}=$ ', num2str(Tign), ' K'], 'Interpreter', 'Latex')
xlabel('x (in m)', 'Interpreter', 'Latex', 'FontSize', 18)
ylabel('time (in s)', 'Interpreter', 'Latex', 'FontSize', 18)


fig4.PaperPosition = [0 0 6 13];
fig4.PaperSize = [6, 13];
orient(fig4,'landscape')
print(fig4,[plot_dir, 'cdf_ign.pdf'],'-dpdf')

