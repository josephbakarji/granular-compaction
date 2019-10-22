function plot_var(Un, plotit, plotfrq, plotmode, xx, t, tend, tsteps, Ti) 
% Plots variables during simulation according to plotmode.
  
global pinf1 pinf2 gam1 gam2

if(plotit == 1)

	if(mod(t, plotfrq) == 0 || isnan(t))

		pps = pk(pinf1, gam1, rho1(Un), e1(Un));
		ppg = pk(pinf2, gam2, rho2(Un), e2(Un));

		if plotmode == 0
			subplot(3, 4, 1); plot(xx, Un(1,:), '.'); title('\alpha_s \rho_s')
			subplot(3, 4, 2); plot(xx, Un(2,:), '.'); title('\alpha_g \rho_g')
			subplot(3, 4, 3); plot(xx, Un(3,:), '.'); title('\rho u')
			subplot(3, 4, 4); plot(xx, Un(4,:), '.'); title('\rho E')
			subplot(3, 4, 5); plot(xx, Un(5,:), '.'); title('\alpha_s')
			subplot(3, 4, 6); plot(xx, Un(6,:), '.'); title('\alpha_s \rho_s e_s')
			subplot(3, 4, 7); plot(xx, Un(7,:), '.'); title('\alpha_g \rho_g e_g')
			subplot(3, 4, 8); plot(xx, Un(8,:), '.'); title('p')
			subplot(3, 4, 9); plot(xx, rho1(Un), '.');	   title('\rho_s');
			subplot(3, 4, 10); plot(xx, e1(Un), '.');	   title('e_s');
			subplot(3, 4, 11); plot(xx, Ti, '.');	   title('Ti');
			subplot(3, 4, 12); plot(xx, (pps - ppg), '.'); title('(ps - pg)');
		elseif plotmode == 1

			subplot(2, 4, 1); plot(xx, a1(Un), '.'); title('\alpha_1')
			subplot(2, 4, 2); plot(xx, pps, '.'); title('p_1')
			subplot(2, 4, 3); plot(xx, ppg, '.'); title('p_2')
			subplot(2, 4, 4); plot(xx, u(Un), '.'); title('u')
			subplot(2, 4, 5); plot(xx, rho1(Un), '.'); title('\rho_1')
			subplot(2, 4, 6); plot(xx, rho2(Un), '.'); title('\rho_2')
			subplot(2, 4, 7); plot(xx, e1(Un), '.'); title('e_1')
			subplot(2, 4, 8); plot(xx, e2(Un), '.'); title('e_2')
		end
	drawnow
	end

end
