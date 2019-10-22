function [FA, vars] = readtxtfile(filename, filedir)
% Reads data from filename.txt from directory filedir and stores results in 3D matrix FA
% FA is with dimensions (size(U, 1), size(U, 2), length(time)), where U is (1+7+compute_temp, length(xx))
	
	%filedir = '../data/'; % FIX FOR MC!!
	temp = strsplit(filename, '.');
	fileprefix = char(temp(1));
	matfile =[filedir, fileprefix, '.mat']; 
	if(~exist(matfile, 'file'))
		disp('mat file does NOT exist')

		fid = fopen([filedir, filename]);
		makevec = @(line) str2num(char(strsplit(line, ' \t')));

		i = 0;
		tline = fgetl(fid);

		vars = makevec(fgetl(fid));
		xend = vars(1);
		ncells = vars(2);
		tend = vars(3);
		tsteps = vars(4);
%		no_relax = vars(5);
%		c1 = vars(6);
%		c2 = vars(7);
%		gam1 = vars(8);
%		gam2 = vars(9);
%		pinf1 = vars(10);
%		pinf2 = vars(11);
%		updatemesh = vars(12);
		if(length(vars)>12)
%			bc = vars(13);
			vbc = vars(14);
%			compute_temp = vars(15);
			temp_mode = vars(16);
%			opt_tol = vars(17);
%			newton_opt = vars(18);
%			cv1 = vars(19);
%			muacc = vars(20);
%			Ls = vars(21);
%			n = vars(22);
%			as0 = vars(23);
%			Ri0 = vars(24);
%			avg_alpha = vars(25);
%			std_alpha = vars(26);
%			distribution = vars(27);
			savemode = vars(28);
		end
%			
%			params = [xend, ncells, tend, tsteps, no_relax, c1, c2, gam1, gam2, pinf1, pinf2, updatemesh, bc, compute_temp, opt_tol, newton_opt, cv1, muacc, Ls, n, as0, Ri0, avg_alpha, std_alpha, distribution];
%		else
%			params = [xend, ncells, tend, tsteps, no_relax, c1, c2, gam1, gam2, pinf1, pinf2, updatemesh];
%		end
%

		varnames = fgetl(fid);
		varnum = length(strsplit(varnames)) - 1; % extra tab (remove in future)
		tline = fgetl(fid);
		j = 1;
		while(ischar(tline))

			A = zeros(varnum, ncells);
			for i = 1:ncells
				tline = fgetl(fid);
				v = makevec(tline);
				try
					A(:,i) = v';
				catch
					keyboard;
				end
			end

			FA(:,:,j) = A;
			j = j + 1;
			if(mod(j, floor(tsteps/10)) == 0)
				disp([num2str(j/tsteps*100), ' % complete'])
			end

			tline = fgetl(fid);
			if(~strcmp(tline, ''))
				disp('done')
				if(~ischar(tline))
					break;
				end
			end

		end

		fclose(fid);
		save([filedir, fileprefix]);
	else
		%disp('mat file exists')
    		load([filedir, fileprefix], 'FA', 'vars');
	end
	    
end

