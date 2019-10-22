function r = E(w)
	if(isvector(w))
		r = w(4)./(w(1)+w(2));		% Total Energy
	else
		r = w(4,:)./(w(1,:) + w(2,:));
	end
end	
