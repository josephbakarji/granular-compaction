function r = u(w)
	if(isvector(w))
		r = w(3)./(w(1) + w(2));	% Velocity	
	else
		r = w(3,:)./(w(1,:) + w(2,:));	% Velocity	
	end
end
