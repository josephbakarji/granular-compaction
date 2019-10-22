function add_extrafiles(foldername, mc_count_des)
% add extra realizations to a simple monte carlo simulation with directory foldername, until mc_count_des files are reached. 

    global gam1 gam2 pinf1 pinf2 c1 c2
    
    addpath('../runfiles/')
    pathfile
    addpath(result_script_dir) % add path of helper functions
    filedir = [MonteCarlo_dir, foldername, '/'];

    if(~exist(filedir, 'file'))
        disp('folder does not exist')
        return
    end

    files = dir([filedir, '*']);

    for i = 1:length(files)
        h = strsplit(files(i).name, '_');
        h{1}
        if(strcmp(h{1}, 'MC'))
            [FA, vars] = readtxtfile(files(i).name, filedir);
            break;
        end
    end

    xend 	= vars(1);
    ncells 	= vars(2);
    tend 	= vars(3);
    tsteps 	= vars(4);
    c1 		= vars(5);
    c2 		= vars(6);
    gam1 	= vars(7);
    gam2 	= vars(8);
    pinf1 	= vars(9);
    pinf2 	= vars(10);
    updatemesh 	= vars(11);
    bc 		= vars(12);
    vbc 	= vars(13);
    compute_temp= vars(14);
    temp_mode 	= vars(15);
    opt_tol 	= vars(16);
    cv1 	= vars(17);
    muacc 	= vars(18);
    Ls 		= vars(19);
    n 		= vars(20);
    as0 	= vars(21);
    Ri0 	= vars(22);
    avg_alpha 	= vars(23);
    std_alpha 	= vars(24);
    distribution= vars(25);
    savemode 	= vars(26);

    savedir 	= filedir;
    plotit 	= 0; 		% plot if 1.
    plotmode 	= 0;
    plotfrq 	= 20;
    printfrq 	= 100;
    varformat 	= '%8.5g';

    % Add to saved data
    u0 			= vbc;        % velocity in m/s
    p0 			= 1e6;       % pressure in pa
    rho1_0 		= 1903;  % hmx density in kg/m^3
    rho2_0 		= 1;     % air density in kg/m^3

    time 		= linspace(0, tend, tsteps);
    xx 			= linspace(0, xend, ncells);
    u0vec 		= zeros(1, length(xx)); u0vec(1) = u0;
    p0vec 		= p0 .* ones(1, length(xx));

    simID = 0;
    simIDmax = 0;
    mc_count_exist = 0;
    for i = 1:length(files)
        basename = strsplit(files(i).name, '.');
        if(strcmp(char(basename(end)), 'txt'))
            temp = strsplit(char(basename{1}), '_');
            simID = str2num(char(temp(end)));
            if(simID >= simIDmax)
                simIDmax = simID;
            end
            mc_count_exist = mc_count_exist + 1;
        end
    end
    mc_count 	= mc_count_des - mc_count_exist;
    nameprefix = 'MC';



    for i = 1:mc_count
        disp(['simulation n.: ', num2str(i)])
        alpha1_0 = get_initial_por(distribution, avg_alpha, std_alpha, xx);
        savefilename = makefilename(distribution, i + simIDmax, nameprefix);

	Uo = make_uo(xx, alpha1_0, u0vec, p0vec, rho1_0, rho2_0, gam1, gam2, pinf1, pinf2);
	solve1D(savefilename, savedir, savemode, compute_temp, temp_mode, plotit, plotfrq, printfrq, plotmode,...
		Uo, time, xx, bc, vbc, updatemesh, varformat, opt_tol, gam1, gam2, pinf1, pinf2, c1, c2, cv1,...
		muacc, Ls, n, as0, Ri0, avg_alpha, std_alpha, distribution)
    end

end
