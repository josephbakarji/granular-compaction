function alpha1_0 = get_initial_por(distribution, avg_alpha, std_alpha, xx, alpha_max)
% Compute initial solid volume fraction based on given probability distirbution
%
%%% Inputs: 
% distribution = 0: Constant value = avg_alpha.
% disttibution = 1: Uniform distirbution with mean avg_alpha, and std std_alpha
% distribution = 2: Gassian distirbution with mean avg_alpha, and std std_alpha
% distribution = 3: Gamma distribution with shape avg_alpha, and scale std_alpha
% distribution = 4: Speific shapes with mean avg_alpha (square, triangle etc.)
% avg_alpha: average porosity
% std_alpha: standard deviation of distirbution.
% xx: 1D coordinates of grid.
%
%%% Output:
% alpha1_0: solid volume fraction of dim length(xx). 
%%%%%%%%%%%%%%%%%%%%%%%%%%

% rng(100)

	if distribution == 0 	%% Deterministic 
		alpha1_0 = avg_alpha .* ones(1, length(xx));

		
	elseif distribution == 1   %% Uniform pore distribution
		for i = 1:length(xx)
			alpha1_0(i) = avg_alpha + std_alpha * (rand-0.5);
		end

	elseif distribution == 2  %% Gaussian pore distribution
		alpha_ciel = 0.98;
		for i = 1:length(xx)
			a = normrnd(avg_alpha, std_alpha);
			if(a >= 1)
				alpha1_0(i) = alpha_ciel;
				disp(['gaussian sampled porosity above limit; setting alpha = ', num2str(alpha_ciel)])
			else
				alpha1_0(i) = a;
			end
		end
	elseif distribution == 3 %% Gamma distribution
		if(nargin == 4)
			alpha_max = 0.98;
		end

		alpha_ciel = 0.98;
		alpha_flr = 0.0001;
		shape_alpha = avg_alpha; % alpha
		scale_alpha = std_alpha; % beta
		for i = 1:length(xx)
			a = alpha_max - gamrnd(shape_alpha, scale_alpha);
			if(a >= alpha_ciel || a <= 0 )
				if( a>= alpha_ciel)
					alpha1_0(i) = alpha_ciel;
					disp(['gamma sampled porosity above limit; setting alpha = ', num2str(alpha_ciel)])
				elseif(a <= 0)
					alpha1_0(i) = alpha_flr;
					disp(['gamma sampled porosity below limit; setting alpha = ', num2str(alpha_flr)])
				end
			else
				alpha1_0(i) = a;
			end
		end
		
	elseif distribution == 4  %% Specific shape

		x1 = xx(end)/5;
		x2 = xx(end)/4;
		x3 = x2 + (x2-x1);
		[~, indx1] = min(abs(xx - x1));
		[~, indx2] = min(abs(xx - x2));
		[~, indx3] = min(abs(xx - x3));
		for i = 1:length(xx)
			if(xx(i)<= xx(indx1))
				alpha1_0(i) = alphamean;
			elseif(xx(i)>xx(indx1) && xx(i) <= xx(indx2))
				alpha1_0(i) = alpha1_0(indx1) + (alphapore - alphamean)/(xx(indx2) - xx(indx1)) * (xx(i) - xx(indx1));
			elseif(xx(i)> xx(indx2) && xx(i) <= xx(indx3))
				alpha1_0(i) = alpha1_0(indx2) + (alphamean - alphapore)/(xx(indx3) - xx(indx2)) * (xx(i) - xx(indx2));
			elseif(xx(i)> xx(indx3))
				alpha1_0(i) = alphamean;
			end
		end

	end

end
