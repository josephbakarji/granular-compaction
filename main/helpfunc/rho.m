function r = rho(w)
	if(isvector(w))	
		r = w(1)+w(2);             % Total Density
	else
		r = w(1, :) + w(2, :);
	end
end
