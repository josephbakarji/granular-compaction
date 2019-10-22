function r = y1(w)
    if(isvector(w))
        r =  w(1)/(w(1)+w(2));
    else
        r = w(1,:)./(w(1,:) + w(2,:));
    end
end
