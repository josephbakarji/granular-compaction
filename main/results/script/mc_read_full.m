function [Tarray, avg_alpha, std_alpha, xx, time, mfrac_s] = mc_read_full(foldername)
% Read Monte Carlo simulation results in foldername and combine them in histmat matrix

pathfile
addpath([maindir, 'helpfunc/']) % add path of helper functions

filedir = [MonteCarlo_dir, foldername, '/'];
files = dir([filedir, '*']);

[FA, params] = readtxtfile(files(end-5).name, filedir);
ncells = params(2);
tend = params(3);
tsteps = params(4);
avg_alpha = params(25);
std_alpha = params(26);
savemode = params(28);

time = linspace(0, tend, tsteps+1); % one added because Initial condition included (?)


% Choose relevant files
relf_num = 0;
for i = 1:length(files)
	h = strsplit(files(i).name, '_');
	fmt = strsplit(h{end}, '.');
	if(strcmp(h{1}, 'MC') && strcmp(fmt{2}, 'txt') && ~strcmp(fmt{1}, '1'))
            relf_num = relf_num + 1;
            rel_files{relf_num} = files(i).name;
        end
end

Tarray = zeros(ncells, tsteps, relf_num); 

disp(['----- Folder name: ', foldername, ' -----'])
percomp_prev = 0;
for i = 1:relf_num
	%disp(files(i).name)

	% Print progress
	percomp = floor(i/length(rel_files)*100);
	if(mod(percomp, 20) == 0)
		if(percomp ~= percomp_prev)
			displayname = [num2str(percomp), ' % complete'];
			disp(displayname);
			percomp_prev = percomp;
		end
	end

	[FA, ~] = readtxtfile(rel_files{i}, filedir);
	for tt = 1:size(FA, 3)
		[A, xx, T] = fa2a2(FA, tt, savemode);
		Tarray(:, tt, i) = T;

		if(savemode == 3)
			mfrac_s(:, tt, i) = A;
		else
			disp('savemode not 3 (fix)')
			mfrac_s(:, tt, i) = nan;	
		end
	end
end
    
sf = strsplit(foldername, '/');
savematname = char(sf(1));
save([MonteCarlo_dir, savematname, '_wlams'])
end

