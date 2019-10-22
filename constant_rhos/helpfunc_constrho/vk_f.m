function r = vk_f(gamk, pinfk, p0k, vk0, p)
	r = vk0 .* (p0k + gamk.*pinfk + (gamk - 1).*p)./(gamk.*(pinfk + p)); % Hugoniot-based relaxation specific volume (Favrie 2013 sec 6b)
end
