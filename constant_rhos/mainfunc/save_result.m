function save_result(savefilename, savedir, savemode, U, xx, varformat, append_flag, compute_temp, temp_mode, T,...
	       	xend, ncells, tend, tsteps, no_relax, c1, c2, gam1, gam2, pinf1, pinf2, updatemesh,...
		bc, vbc, opt_tol, newton_opt, cv1, muacc, Ls, n, as0, Ri0, avg_alpha, std_alpha, distribution)

% Saves results to savefilename.txt in savedir folder according to savemode format.

	% TODO: REMOVE varset_save?
	savefilename_full = [savedir, savefilename];

	if append_flag == 0 

		keyset = {'xend', 'ncells', 'tend', 'tsteps', 'no_relax', 'c1', 'c2', 'gam1', 'gam2', 'pinf1', 'pinf2', 'updatemesh', 'bc', 'v_bc', 'compute_temp', 'temp_mode', 'opt_tol', 'newton_opt', 'cv1', 'muacc', 'Ls', 'n', 'as0', 'Ri0', 'avg_alpha', 'std_alpha', 'distribution', 'savemode'};
		keytypes={ '%5f', '%d', '%g', '%d', '%d', '%d', '%d', '%4.2f', '%4.2f', '%g', '%2.1f', '%d', '%d', '%g' ,'%d', '%d', '%g', '%d', '%5.2f', '%4g', '%f', '%d', '%f', '%g', '%4.2f', '%4.2f', '%d', '%d'};
		valueset = [xend, ncells, tend, tsteps, no_relax, c1, c2, gam1, gam2, pinf1, pinf2, updatemesh, bc, vbc, compute_temp, temp_mode, opt_tol, newton_opt, cv1, muacc, Ls, n, as0, Ri0, avg_alpha, std_alpha, distribution, savemode];
		optstr = '';
		typestr = ''; 
		varstr1 = '';
		vartype1 = '';

		for i  = 1:length(valueset)
			optstr = strcat(optstr, keyset{i}, ' \t ');
			typestr = strcat(typestr, keytypes{i}, ' \t ');
		end

		[As, varset_save] = make_array(xx, U, T, savemode);

		for i = 1:length(varset_save)
			varstr1 = strcat(varstr1, varset_save{i}, ' \t ');
			vartype1 = [vartype1, varformat, '\t'];
		end

		optstr = strcat(optstr, ' \n');
		typestr = strcat(typestr, ' \n');
		varstr1 = strcat(varstr1, ' \n\n');
		vartype1 = strcat(vartype1, '\n');


		if ~isempty(savefilename)
			fileid = fopen(savefilename_full, 'w');
			fprintf(fileid, optstr);
			fprintf(fileid, typestr, valueset);
			fprintf(fileid, varstr1);
			fprintf(fileid, vartype1, As);
			fclose(fileid);
		end
	else

		[As, varset_save] = make_array(xx, U, T, savemode);
		varnum = length(varset_save);	% FIX: make function of varidx

		vartype1 = '';
		for i  = 1:varnum
			vartype1 = [vartype1, varformat, '\t'];
		end
		vartype1 = strcat(vartype1, '\n');

		if ~isempty(savefilename)
			[fileid, errmsg] = fopen(savefilename_full, 'a');
			if(~isempty(errmsg))
				disp(errmsg)
				keyboard	
			end
			while(~exist('fileid')), end
			fprintf(fileid, '\n');
			fprintf(fileid, vartype1, As);
			fclose(fileid);
		end
	end
end
