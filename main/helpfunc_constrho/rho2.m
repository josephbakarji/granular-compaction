function r = rho2(w)
	if(isvector(w))
		r  = w(2)./(1-w(5));		% Density phase 2
	else
		r = w(2,:) ./(1-w(5,:));
	end
end
