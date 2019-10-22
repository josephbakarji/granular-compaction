function r = rho1(w)
	if(isvector(w))
		r = w(1)./w(5);            % Density phase 1
	else
		r = w(1,:)./w(5,:);            % Density phase 1
	end
end
