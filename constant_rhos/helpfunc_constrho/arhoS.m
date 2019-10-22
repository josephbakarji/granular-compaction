function r = arhoS(arho, Sk, Sm, uk)
	r =  arho * (Sk - uk)/(Sk - Sm);
end

