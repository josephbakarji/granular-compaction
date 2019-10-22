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

    xend = vars(1);
    ncells = vars(2);
    tend = vars(3);
    tsteps = vars(4);
    no_relax = vars(5);
    c1 = vars(6);
    c2 = vars(7);
    gam1 = vars(8);
    gam2 = vars(9);
    pinf1 = vars(10);
    pinf2 = vars(11);
    updatemesh = vars(12);
    bc = vars(13);
    vbc = vars(14);
    compute_temp = vars(15);
    temp_mode = vars(16);
    opt_tol = vars(17);
    newton_opt = vars(18);
    cv1 = vars(19);
    muacc = vars(20);
    Ls = vars(21);
    n = vars(22);
    as0 = vars(23);
    Ri0 = vars(24);
    avg_alpha = vars(25);
    std_alpha = vars(26);
    distribution = vars(27);
    savemode = 3;

    etam            = 2e-3; % Missing in saved data (not used)
    Y1          	= 400e6; 
    B           	= 5765;
    Tm          	= 1000; % ??? random

    savedir 	= filedir;

    plotit 		= 0; 		% plot if 1.
    plotmode 	= 0;
    plotfrq 	= 20;
    printfrq 	= 100;
    varformat = '%8.5g';

    u0 		= 60;        % velocity in m/s
    p0 		= 1e6;       % pressure in pa
    rho1_0 		= 1903;  % hmx density in kg/m^3
    rho2_0 		= 1;     % air density in kg/m^3

    time 		= linspace(0, tend, tsteps);
    xx 		= linspace(0, xend, ncells);
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


    %% problem parameters


    for i = 1:mc_count

        disp(['simulation n.: ', num2str(i)])
        alpha1_0 = get_initial_por(distribution, avg_alpha, std_alpha, xx);

        savefilename = makefilename(distribution, i + simIDmax, nameprefix);

        Uo = make_uo(alpha1_0, p0vec, u0vec, rho1_0, rho2_0, gam1, gam2, pinf1, pinf2);
        saur_main(savefilename, savedir, savemode, no_relax, compute_temp, temp_mode, plotit, plotfrq, printfrq, plotmode,...
            Uo, time, xx, bc, u0, updatemesh, varformat, opt_tol, newton_opt,...
            gam1, gam2, pinf1, pinf2, c1, c2, cv1, muacc, Ls, n, as0, Ri0, etam, Y1, B, Tm,...
            avg_alpha, std_alpha, distribution);


    end

end
