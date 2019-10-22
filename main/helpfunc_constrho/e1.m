function r = e1(w)
	if(isvector(w))
		r =  w(6)./w(1);            % Internal energy phase 1
	else
		r =  w(6,:)./w(1,:);            % Internal energy phase 1
	end
end
