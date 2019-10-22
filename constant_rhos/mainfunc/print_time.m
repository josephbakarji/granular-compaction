function print_time(t, tsteps, tend, printfrq)
	if(mod(t, printfrq)==0)
		str= sprintf('time = %4.2e / %4.2e ms \t -\t time step: %d / %d',  t/tsteps*tend*1000, tend*1000, t, tsteps);
		disp(str);
	end
end

