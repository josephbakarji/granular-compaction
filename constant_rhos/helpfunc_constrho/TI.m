function [r, Tps] = TI(Tc, dt, Rint, betas, alphag, n, mu, Ls, ps, pg)
	Tps = (dt - Rint) .* mu .* (Rint .* betas)./ (3 .* alphag .* n .* Ls) .* (ps - pg);
	r =  Tc + Tps;
end
