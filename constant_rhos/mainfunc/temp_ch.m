function Ti = temp_ch(Uo, Un, dt, pn, Ri, N2, cv1, Y1, etam, B, Tio, Tm)
% Compute temperature dissipation based on CH model (see Nesterenko 2013 "Dynamics of Heterogeneous Materials" p. 315)


Y = @(T) heaviside(Tm - T) .* Y1 .* (1 - T./Tm);
eta = @(T) etam .* exp(B .* (1./T + 1/Tm));


for i = 1:size(Un, 2)
	dadt = 1/dt .* (3/(4*pi * N2)).^(1/3) .* (a2(Un(:, i))^(1/3) - a2(Uo(:, i))^(1/3));
	Ti_pdiss = Tio(i) + dt .* ( - 2/(rho1(Un(:,i)) .* cv1) .* Y(Tio(i)) .* dadt ./ a2(Un(:,i)) + 12 * eta(Tio(i)) .* dadt.^2 ./ a2(Un(:,i)).^2);
	Ti_diff = 0;
	Ti_react = 0;
	Ti(i) = Ti_pdiss + Ti_diff;
end





end



