function [Ti, Rim1] = temperature_simple(Un, Rim1, dt, es0, ps0, vs0, N2, cv1, muacc, Ls, as0, c1, c2)

	% Compute temperature according to known diffusion profile in spherical shell.
	% See Saurel 2018 "Multiscale multiphase modeling of ..." for more details.
	% TODO:
	% - Fix issue with dt0 sometimes giving Nan, or out of range for n = 1.
	% - Try for n = 5.

	global dt0f_count

	Tav = @(P) 0.4385 .* (P.*1e-9).^2 + 15.18 .* (P.*1e-9) + 310; % ORIGINAL P in Pa; from appendix

	% MAIN Variables
	Ri = ( 3 .* a2(Un) ./ (N2 * 4 * pi) ).^(1/3);
	Re = Ri ./ a2(Un).^(1/3); % for deterministic Ri
	Tav1 = Tav( Un(8,:) );
	Tc = 1/cv1 .* (es0 - 0.5 .* (ps0 + Un(8,:)) .* (1./(rho1(Un) - vs0) ) ); % REMOVE 200! Verify!

	% Temporary Fix
	Z1 = rho1(Un) .* c1;
	Z2 = rho2(Un) .* c2;
	mu = muacc .* (4 .* pi .* Ri.^2 * N2./(Z1 + Z2));

	for i = 1:size(Un, 2)
		dBvec(i) = dBconf(Un(5,i), as0);
	end
	betas = Un(1,:) .* dBvec;

	dt0f_c = 1;
	flag = 0;
	% Compute approx. interface temperature (using n = 1)
	for j = 1:size(Un, 2)
		dtopt = @(dtbl) optexp_new(Tc(i), Tav1(i), dtbl, dt, Ri(i), Rim1(i), Re(i), betas(i), 1, Ls);
		if(j > 1)
			try
				dtbl(j) = fzero(dtopt, [0, dtbl(j-1) + 1]); % MAKE MORE ROBUST!
			catch
				if(flag ~= 0 )
					flag = 1;
					disp('negative dtbl at j > 1')
					dt0f_count = dt0f_count + 1;
				end   
				dt0f_c = dt0f_c + 1;
				dtbl(j) = (Re(j) + Ri(j))/2 + Ri(j);

			end
		else
			try
				dtbl(j) = fzero(dtopt, [0, (Re(j) + Ri(j))/2 + 1]);
			catch
				dtbl(j) = (Re(j) + Ri(j))/2 + Ri(j);
				if(flag ~= 0 )
					flag = 1;
					disp('negative dt0 at j = 1')
					dt0f_count = dt0f_count + 1;
				end    
			end	
		end

	end
	if(flag == 1)
		disp(['n. dt0 TOTAL fault = ', num2str(dt0f_count), ...
			' - n. dt0 LOCAL fault = ', num2str(dt0f_c), ...
		       	' -  With Un size ', num2str(size(Un, 2))])
	end
	
	%keyboard
	[Ti, Tps] = TI2(Tc, dtbl, Ri, Rim1, betas, 1, Ls, dt);
	Ti(1) = Ti(2);
		
	disp(['Tc = ', num2str(Tc(1:3))]);
	disp(['Tps = ', num2str(Tps(1:3))]);
	disp(['Ti = ', num2str(Ti(1:3))]); 
	disp(['Rim1 = ', num2str(Rim1(1:3))]); 
	disp(['Ri = ', num2str(Ri(1:3))]); 
	
	Rim1 = Ri;




	% Use n = 5
	%	for i = 1:size(Un, 2)
	%		if dt0(i)<(Re(i) - Ri(i)) && dt0(i)>0
	%		    Ti(i) = Tc(i) + (dt0(i) - Ri(i)) .* mu(i) .* Ri(i) .* Un(1,i).*dBconf(Un(5,i), as0) ./ (3 .* (1 - Un(5,i)) .* n .* Ls) .* (pps(i) - ppg(i)); %dBconf0
	%		elseif dt0(i)>(Re(i) - Ri(i))
	%		    Y = - Re(i).*Ri(i) - (Re(i) - Ri(i)).^2./3;
	%		    Z = - Re(i).^2/(n+1) + 2*Re(i)*(Re(i) - Ri(i))./(n+2) - (Re(i) - Ri(i)).^2/(n+3);
	%		    Ti(i) = ( Ls*n.*Tav1(i).*(Re(i).^3 - Ri(i).^3)./(3.*(Ri(i) - Re(i)).*(Y-Z)) +...
	%			 (Re(i) - Ri(i)) .* mu(i) .* Ri(i) .* Un(1,i).*dBconf(Un(5,i), as0) ./ (3 .* (1 - Un(5,i)) .* n .* Ls) .* (pps(i) - ppg(i)) )...
	%			./( (n * Ls .* Y)./(Y - Z));
	%		elseif dt0(i)<0
	%		    Ti(i) = Tc(i) + (dt0(i) - Ri(i)) .* mu(i) .* Ri(i) .* Un(1,i).*dBconf(Un(5,i), as0) ./ (3 .* (1 - Un(5,i)) .* n .* Ls) .* (pps(i) - ppg(i)); %dBconf0
	%		    disp('Radius too small?')
	%		end
	%	end

end 
