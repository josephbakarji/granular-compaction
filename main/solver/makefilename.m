function filename = makefilename(CASE, simID, casename)
	if nargin == 2
		filename =sprintf(['case_%d_simID_%d.txt'], CASE, simID);
	elseif nargin == 3 
		filename =sprintf([casename,'_%d_simID_%d.txt'], CASE, simID);
	end

end
