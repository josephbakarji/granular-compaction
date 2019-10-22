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
		if(length(vars)>12)
			vbc = vars(13);
			temp_mode = vars(15);
			savemode = vars(26);
		end
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

