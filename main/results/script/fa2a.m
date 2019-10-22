function [A, xx, T] = fa2a(FA, tt, compute_temp)
% Transform the full space-time tensor FA to A at time tt.
% If compute_temp == 1, last row is temperature.

	if(compute_temp)
		xx = FA(1,:,tt);
		A = FA(2:end-1,:,tt);
		T = FA(end, :, tt);
	else
		xx = FA(1,:,tt);
		A = FA(2:end,:,tt);
	end
end
