function r = pI(rho1, c1, p1, rho2, c2, p2)
	r = (rho2.* c2 .* p1 + rho1 .* c1 .* p2) ./ (rho1 .* c1 + rho2 .* c2);
end
