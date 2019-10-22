function r = pk(pinfk, gamk, rhok, ek) 
	r = rhok .* ek .* (gamk - 1) - gamk .* pinfk; 	% Stiffened gas phase pressure
end
