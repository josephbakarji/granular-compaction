function r = a1(w)
	if(isvector(w))
		r = w(5);                  % Volume fraction phase 1
	else
		r = w(5,:);
	end
end
