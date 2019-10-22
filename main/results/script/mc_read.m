function [histmat, Tdis, xx, time] = mc_read(foldername)
% Read Monte Carlo simulation results in foldername and combine them in histmat matrix

pathfile
addpath([maindir, 'helpfunc/']) % add path of helper functions
filedir = [MonteCarlo_dir, foldername, '/'];
histnum = 200;
TmaxMargin = 450;
TminMargin = 5;
files = dir([filedir, '*']);

[FA, params] = readtxtfile(files(end-5).name, filedir);
ncells = params(2);
tend = params(3);
tsteps = params(4);
time = linspace(0, tend, tsteps+1); % one added because Initial condition included (?)
[A, xx, T] = fa2a2(FA, 50, 2);


Tmin = 200; %min(T) - TminMargin;
Tmax = 1400; %max(T) + TmaxMargin;

Tdis = linspace(Tmin, Tmax, histnum);
histmat = zeros(size(FA, 2), size(FA, 3), histnum);

disp(['----- Folder name: ', foldername, ' -----'])
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
            
            [FA, params] = readtxtfile(files(i).name, filedir);
            for tt = 1:size(FA, 3)
                [A, xx, T] = fa2a2(FA, tt, 2);
                for j = 1:length(xx)
                    Tgrid = find(Tdis < T(j), 1, 'last');
                    histmat(j, tt, Tgrid) = histmat(j, tt, Tgrid) + 1;
                end
            end
	end
end
clear FA params
sf = strsplit(foldername, '/');
savematname = char(sf(1));
save([MonteCarlo_dir, savematname])

end

