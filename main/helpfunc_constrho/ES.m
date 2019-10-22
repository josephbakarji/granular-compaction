function r = ES(rhok, Ek, uk, Sk, pk, rhoSk, Sm, pS)
	r = ( rhok * Ek * (uk - Sk) + pk*uk - pS*Sm)/(rhoSk*(Sm - Sk));
end
