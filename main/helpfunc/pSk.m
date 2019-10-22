function r = pSk(gamk, pk, pinfk, rhoSk, rhok) 
	r = (pk + pinfk) * ( (gamk - 1)*rhok - (gamk + 1)*rhoSk )/( (gamk - 1)*rhoSk - (gamk + 1)*rhok ) - pinfk;
end
