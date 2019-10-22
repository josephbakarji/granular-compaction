function r = rho1(w)
	global rhos
	if(isvector(w))
		r = rhos;  % Density phase 1
	else
		r = rhos .* ones(size(w(1, :))); % Density phase 1
	end
end
