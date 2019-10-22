function r = y2(w)
	if(isvector(w))
        r =  w(2)/(w(1)+w(2));
    else
        r = w(2,:)./(w(1,:) + w(2,:));
    end
end
