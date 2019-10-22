function [As, varset_save] = make_array(xx, U, T, savemode)
% builds array to be saves in text file according depending on savemode.

% TODO:
% - add possibility to select varset_save as input (choose variables)
% - add possibility to save pure variables: xx, a1, rho1, rho2, u, E etc.
	if(savemode == 0)
		As = [xx; U];
		varset_save = {'xx', 'a1*rho1', 'a2*rho2', 'rho*u', 'rho*E', 'a1', 'a1*rho1*e1', 'a2*rho2*e2', 'p'};
	elseif(savemode == 1)
		As = [xx; U; T];
		varset_save = {'xx', 'a1*rho1', 'a2*rho2', 'rho*u', 'rho*E', 'a1', 'a1*rho1*e1', 'a2*rho2*e2', 'p', 'T'};
	elseif(savemode == 2)
		As = [xx; T];
		varset_save = {'xx', 'T'};
	elseif(savemode == 3) % Added to account for solid mass fraction in reaction rate (Bdzil)
		lam_s = y1(U);
		As = [xx; lam_s; T];
		varset_save = {'xx', 'lam_s', 'T'};
	end
end
