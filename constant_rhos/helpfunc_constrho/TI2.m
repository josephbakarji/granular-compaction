function [r, Tps] = TI2(Tc, dtbl, Rint, Rim1, betas, n, Ls, dt)
	Tps = - betas .* (Rint - Rim1)./dt .* (dtbl - Rint) ./ (n .* Ls);
	r =  Tc + Tps;
end
