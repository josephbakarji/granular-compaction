function r = a2(w)
	if(isvector(w))
		r = 1 - w(5);                % Volume fraction phase 2
	else
		r = 1 - w(5,:);
	end
end

