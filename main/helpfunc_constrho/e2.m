function r = e2(w)
	if(isvector(w))
		r =  w(7)./w(2);            % Internal energy phase 1
	else
		r =  w(7,:)./w(2,:);            % Internal energy phase 1
	end
end
