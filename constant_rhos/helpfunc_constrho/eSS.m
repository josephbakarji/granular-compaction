function r = eSS(gam, pinf, pSkk, rhoSk)
	r = (pSkk + gam * pinf)/((gam - 1)*rhoSk);
end
