function [Tsup_time, xx, time, avg_alpha, std_alpha] = mc_read_timecorr(foldername, Tig)
	% Read Monte Carlo simulation results in foldername and combine them in histmat matrix

	pathfile
	addpath([maindir, 'helpfunc/']) % add path of helper functions
	filedir = [MonteCarlo_dir, foldername, '/'];
	
    



	files = dir([filedir, '*']);

	[FA, params] = readtxtfile(files(end-5).name, filedir);
	tend = params(3);
	tsteps = params(4);
    avg_alpha = params(25);
    std_alpha = params(26);
    savemode = params(28);
	
    dt = tend/tsteps;
	time = linspace(0, tend, tsteps+1); % one added because Initial condition included (?)
	[A, xx, T] = fa2a2(FA, 50, savemode);


	Tsup_time = [];
	percomp_prev = 0;
	for i = 1:length(files)
		h = strsplit(files(i).name, '_');
		fmt = strsplit(h{end}, '.');
		if(strcmp(h{1}, 'MC') && strcmp(fmt{2}, 'txt') && ~strcmp(fmt{1}, '1'))
			%disp(files(i).name)
			percomp = floor(i/length(files)*100);
			if(mod(percomp, 20) == 0)
				if(percomp ~= percomp_prev)
					displayname = [num2str(percomp), ' % complete'];
					disp(displayname);
					percomp_prev = percomp;
				end
			end

			[FA, ~] = readtxtfile(files(i).name, filedir);
			Tmat = [];
			for tt = 1:size(FA, 3)
				[A, xx, T] = fa2a2(FA, tt, savemode);
				Tmat = [Tmat, T'];
			end
			%keyboard	
			Tsi = Tmat > Tig;
			Tsup_time = [Tsup_time, sum(Tsi, 2) * dt]; % APProximation: only accurate if temperatures have a single maximum.
		end
	end



end

