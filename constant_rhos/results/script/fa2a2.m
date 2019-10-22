function [A, xx, T] = fa2a2(FA, tt, savemode)
%  FA to A at tt, adapted to savemode environment.
% TODO:
% Unify fa2a2 with fa2a.

	if(savemode == 1)
		xx = FA(1,:,tt);
		A = FA(2:end-1,:,tt);
		T = FA(end, :, tt);
	elseif(savemode == 0)
		xx = FA(1,:,tt);
		A = FA(2:end,:,tt);
	elseif(savemode == 2)
		A = nan;
		xx = FA(1,:,tt);
		T = FA(2,:,tt);
	elseif(savemode == 3)
		xx = FA(1,:,tt);
		A = FA(2,:,tt); % lam_s
		T = FA(3,:,tt);
	end
end
