function db = dBconf(as, a0)
% Derivative of confugurational energy dB/da, with B(a) is a known function of volume fraction.
% Inputs:
% as: volume fraction
% a0: volume fraction at zero confinement pressure (underformed grains)
% Outputs:
% db: derivative dB/da.

% in Favrie et. al (2013)
m = 1.01;
np = 1.4;
%a = 850*np; 
%a = 60000; % From Saurel for EM (2018)
a = 2000; % Tweaked

% in Saurel et. al (2018)
% np = 1.1;
% a0 = 0.982999;
%a = 1000;


g = @(as) (1-as).^m * log(1-as);
dg = @(as) - m.*(1-as).^(m-1).*log(1-as) - (1-as).^(m-1);

if(as<=1 && as>a0)
    db = a .* (g(as) - g(a0) - dg(a0).*(as-a0)).^(np-1) .* (dg(as) - dg(a0));
elseif(as>=0 && as<=a0)
    db = 0;
else
%     as = 0.999;
%     db = a .* (g(as)-g(a0)-dg(a0).*(as-a0)).^(np-1) .* (dg(as) - dg(a0));

    db = nan;
    disp('\alpha_s not in range')
end
    
