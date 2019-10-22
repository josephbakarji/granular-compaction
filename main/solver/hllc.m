function flux = hllc(Ul, Ur, Uls, Urs, Fl, Fr, Sl, Sr, Sm)
% Compute wave speed based on HLLC algorithm
% For more details see Riemann Solvers and Numerical Methods for Fluid Dynamics, Toro p. 315

    if(Sl >= 0)
        flux = Fl;
    elseif(Sm >= 0 && Sl < 0)
        flux = Fl + Sl * (Uls - Ul);
    elseif(Sr >= 0 && Sm < 0)
        flux = Fr + Sr * (Urs - Ur);
    elseif(Sr <= 0)
        flux = Fr;
    else
        flux = nan;
        disp('none of the hllc conditions met')
    end

end
