clear variables
close all
clc

% File to compare plots manually selected from folder/file below

pathfile
addpath([maindir, 'helpfunc/']) % add path of helper functions

plotm = 0;
plotcompare = 0;
plotcompare2 = 1;
plotsingle = 0;

if plotcompare
	file1 = 'case_3_simID_16.txt';
	file2 = 'case_3_simID_17.txt';
	f1 = 'compare_saur_corr_finegrid_1_simID_1.txt';
	f2 = 'compare_saur_corr_finegrid_2_simID_3.txt';
	f3 = 'compare_saur_corr_coarsegrid_2_simID_4.txt';
	f4 = 'compare_saur_corr_coarsegrid_1_simID_2.txt';
	fp1 = 'compare_saur_pdiff_coarsegrid_1_simID_1.txt';
	fp2 = 'compare_saur_pdiff_coarsegrid_2_simID_2.txt';
	fpn1 = 'compare_saur_pdiffn_coarsegrid_1_simID_1.txt';
	fpn2 = 'compare_saur_pdiffn_coarsegrid_2_simID_2.txt';
	ffpn2 = 'compare_saur_pdiffn_finegrid_2_simID_1.txt';


	pertend = 0.4; % percentage of simulation


	[FA1, params1] = readtxtfile(ffpn2, data_dir);
	[FA2, params2] = readtxtfile(fpn2, data_dir);


	tsteps1 = params1(4);
	tsteps2 = params2(4);
	compute_temp1 = 0;
	compute_temp2 = params2(14);
	tt1 = floor(pertend * tsteps1);
	tt2 = floor(pertend * tsteps2);

	[A1, xx1] = fa2a(FA1, tt1, compute_temp1);
	[A2, xx2] = fa2a(FA2, tt2, compute_temp2);

	% assuming same for both
	gam1 = params1(8);
	gam2 = params1(9);
	pinf1 = params1(10);
	pinf2 = params1(11);

	p1_1 = pk(pinf1, gam1, rho1(A1), e1(A1));
	p2_1 = pk(pinf2, gam2, rho2(A1), e2(A1));
	p1_2 = pk(pinf1, gam1, rho1(A2), e1(A2));
	p2_2 = pk(pinf2, gam2, rho2(A2), e2(A2));

	fig1 = figure;	
        subplot(2, 4, 1); plot(xx1, a1(A1), '.'); hold on;   plot(xx2, a1(A2), '.');    title('\alpha_1')
        subplot(2, 4, 2); plot(xx1, p1_1, '.');	  hold on;   plot(xx2, p1_2, '.');      title('p_1')
        subplot(2, 4, 3); plot(xx1, p2_1, '.');	 hold on;    plot(xx2, p2_2, '.');      title('p_2')
        subplot(2, 4, 4); plot(xx1, u(A1), '.'); hold on;    plot(xx2, u(A2), '.'); title('u')
        subplot(2, 4, 5); plot(xx1, rho1(A1), '.');  hold on;plot(xx2, rho1(A2), '.'); 	title('\rho_1')
        subplot(2, 4, 6); plot(xx1, rho2(A1), '.'); hold on; plot(xx2, rho2(A2), '.'); 	title('\rho_2')
        subplot(2, 4, 7); plot(xx1, e1(A1), '.');  hold on;  plot(xx2, e1(A2), '.');    title('e_1')
        subplot(2, 4, 8); plot(xx1, e2(A1), '.');  hold on;  plot(xx2, e2(A2), '.');    title('e_2')
	legend({'fine', 'coarse'}, 'Location', 'southwest')

	fig1.PaperPositionMode = 'manual';
	orient(fig1,'landscape')
	print(fig1,[plot_dir, 'Fig_pcomp_1.pdf'],'-dpdf')
end




if plotsingle
	
	pertend = 0.95; % percentage of simulation


	filename = 'MC_0_simID_1.txt';
	[FA1, params1] = readtxtfile(filename, data_dir);

	tottsteps = size(FA1, 3);

	tsteps1 = params1(4);
	compute_temp1 = params1(14);
	tt1 = floor(pertend * tottsteps);

	[A1, xx1] = fa2a(FA1, tt1, compute_temp1);

	% assuming same for both
	gam1 = params1(8);
	gam2 = params1(9);
	pinf1 = params1(10);
	pinf2 = params1(11);

	p1_1 = pk(pinf1, gam1, rho1(A1), e1(A1));
	p2_1 = pk(pinf2, gam2, rho2(A1), e2(A1));

	fig1 = figure;	
        subplot(2, 4, 1); plot(xx1, a1(A1), '.');  title('\alpha_1')
        subplot(2, 4, 2); plot(xx1, p1_1, '.');	   title('p_1')
        subplot(2, 4, 3); plot(xx1, p2_1, '.');	   title('p_2')
        subplot(2, 4, 4); plot(xx1, u(A1), '.');   title('u')
        subplot(2, 4, 5); plot(xx1, rho1(A1), '.');title('\rho_1')
        subplot(2, 4, 6); plot(xx1, rho2(A1), '.');title('\rho_2')
        subplot(2, 4, 7); plot(xx1, e1(A1), '.');  title('e_1')
        subplot(2, 4, 8); plot(xx1, e2(A1), '.');  title('e_2')

	fig1.PaperPositionMode = 'manual';
	orient(fig1,'landscape')
	print(fig1,[plot_dir, 'Fig_compaction.pdf'],'-dpdf')
end





if plotcompare2

	pertend = 0.6; % percentage of simulation
	
	f1 = 'singles/SINGLE_2_simID_1.txt';
	f2 = 'singles/SINGLE022_2_simID_1.txt';

	[FA1, params1] = readtxtfile(f1, data_dir);
	[FA2, params2] = readtxtfile(f2, data_dir);


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

	fig1 = figure;
	subplot(2, 2, 1);
	w1 = plot(xx1c, a11(idx0:end), '.'); hold on;
	w2 = plot(xx2c, a12(idx0:end), '.');
	t1 = title('Solid Volume Fraction');
	x1 = xlabel('x (in m)');
	y1 = ylabel('$\alpha_s$');
	legend({'\sigma^2 = 0.02', '\sigma^2 = 0.04'}, 'Location', 'southwest', 'Fontsize', 10)

        subplot(2, 2, 2);
	w3 = plot(xx1c, T1(idx0:end), '.'); hold on;
	w4 = plot(xx2c, T2(idx0:end), '.');
	t2 = title('Pore Surface Temperature ');
	x2 = xlabel('x (in m)');
	y2 = ylabel('$T_I$ (in K)');
	legend({'\sigma^2 = 0.02', '\sigma^2 = 0.04'}, 'Location', 'northeast', 'Fontsize', 12)

        subplot(2, 2, 3);
	w5 = plot(xx1c, p1_1(idx0:end), '.'); hold on;
	w6 = plot(xx2c, p1_2(idx0:end), '.');
	t3 = title('Solid Phase Pressure ');
	x3 = xlabel('x (in m)');
	y3 = ylabel('$p_s$ (in Pa)');
	legend({'\sigma^2 = 0.02', '\sigma^2 = 0.04'}, 'Location', 'northeast', 'Fontsize', 12)

        subplot(2, 2, 4);
	w7 = plot(xx1c, u1(idx0:end), '.'); hold on;
	w8 = plot(xx2c, u2(idx0:end), '.');
	t4 = title('Velocity');
	x4 = xlabel('x (in m)');
	y4 = ylabel('$u$ (in m/s)');
	legend({'\sigma^2 = 0.02', '\sigma^2 = 0.04'}, 'Location', 'northeast', 'Fontsize', 12)


	set([w1, w2, w3, w4, w5, w6, w7, w8], 'MarkerSize', 10)
	set([x1, x2, x3, x4, y1, y2, y3, y4 , t1, t2, t3, t4], 'Interpreter', 'Latex', 'Fontsize', 16)

	fig1.PaperPositionMode = 'manual';
	orient(fig1,'landscape')
	print(fig1,[plot_dir, 'basic_result.pdf'],'-dpdf')

	fig2 = figure;
	subplot(1, 2, 1);
	w1 = plot(xx1c, a11(idx0:end), '.'); hold on;
	w2 = plot(xx2c, a12(idx0:end), '.');
	t1 = title('Solid Volume Fraction');
	x1 = xlabel('x (in m)');
	y1 = ylabel('$\alpha_s$');
	legend({'$\sigma_\xi^2 = 0.02$', '$\sigma_\xi^2 = 0.04$'}, 'Location', 'southwest', 'Interpreter', 'Latex','Fontsize', 16)

        subplot(1, 2, 2);
	w3 = plot(xx1c, T1(idx0:end), '.'); hold on;
	w4 = plot(xx2c, T2(idx0:end), '.');
	t2 = title('Pore Surface Temperature ');
	x2 = xlabel('x (in m)');
	y2 = ylabel('$T_{ps}$ (in K)');
	legend({'$\sigma_\xi^2 = 0.02$', '$\sigma_\xi^2 = 0.04$'}, 'Location', 'northeast', 'Interpreter', 'Latex', 'Fontsize', 16)


	set([w1, w2, w3, w4, w5, w6, w7, w8], 'MarkerSize', 10)
	set([x1, x2, x3, x4, y1, y2, y3, y4 , t1, t2, t3, t4], 'Interpreter', 'Latex', 'Fontsize', 16)

	fig2.PaperPosition = [0 0 5 13];
	fig2.PaperSize = [5, 13];
	orient(fig2,'landscape')
	print(fig2,[plot_dir, 'alpha_temp.pdf'],'-dpdf')


end

