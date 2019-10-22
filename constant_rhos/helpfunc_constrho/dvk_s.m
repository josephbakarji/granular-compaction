function r = dvk_s(gamk, pinfk, pk0, vk0, pk, PI) 
	r = - vk0 .* (pk0 + gamk.*pinfk + (gamk - 1).*PI)./(pk + gamk.*pinfk  + (gamk - 1).*PI).^2;
end
