function Uo = make_uo_fromTemp(xx, alpha1_0, u0vec, p0vec, rho1_0, rho2_0, gam1, gam2, pinf1, pinf2, cv1, cv2, T1_0, T2_0)

	%%% CONSERVATIVE VARIABLES And Fluxes (set in Vectors Uo and Un)
	% w = [    w(1)     ,    w(2)     ,   w(3)      ,  w(4)       ];
	% U = [ (a1*rho1)   , (a2*rho2)   , (rho*u)     , (rho*E)     ];
	% F = [ (a1*rho1)*u , (a2*rho2)*u , (rho*u^2+p) , (rho*E+p)*u ];
	%%%
	%%% NON-CONSERVATIVE VARIABLES
	% w  = [ w(5),  w(6)     ,   w(7)     ,  w(8)]
	% Uv = [ a1 , a1*rho1*e1 , a2*rho2*e2 ,   p  ]
	%
	% Note: w(8) is not computed in main loop.

	Uo = zeros(8, length(xx));    % old value of state
	alpha2_0 = 1 - alpha1_0;
	rho0 = (alpha1_0 .* rho1_0 + alpha2_0 .* rho2_0);
	
	onesvec = ones(1, length(xx));
	e1_0 = T1_0 .* cv1 .* onesvec;
	e2_0 = T2_0 .* cv2 .* onesvec;
	p1_0 = (gam1 - 1).* rho1_0 .* e1_0 - pinf1 .* gam1;
	p2_0 = (gam2 - 1).* rho2_0 .* e2_0 - pinf2 .* gam2;
	p0 = alpha1_0 .* p1_0 + alpha2_0 .* p2_0;
	rhoE0 = alpha1_0 .* rho1_0 .* e1_0 + alpha2_0 .* rho2_0 .* e2_0 + 0.5 .* rho0 .* u0vec.^2; % Fixed 0.5 * rho0 * u^2

	Uo(1,:) = alpha1_0 .* rho1_0;
	Uo(2,:) = alpha2_0 .* rho2_0;
	Uo(3,:) = rho0 .* u0vec;
	Uo(4,:) = rhoE0;
	Uo(5,:) = alpha1_0;
	Uo(6,:) = alpha1_0 .* rho1_0 .* e1_0;
	Uo(7,:) = alpha2_0 .* rho2_0 .* e2_0;
	Uo(8,:) = p0;
	keyboard
end
