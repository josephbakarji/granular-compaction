function r = vk_s(gamk, pinfk, p0, vk0, p, PI)
	r = vk0 .* (p0 + gamk.*pinfk + (gamk - 1).*PI)./(p + gamk.*pinfk  + (gamk - 1).*PI);% (Saurel 2009 III.4-)
end
