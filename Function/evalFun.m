function out = evalFun(mode, n, psf, SBR, Polar, A, phi)
%EVALFUN calculates the loglikelihood function
%   out = evalFun(mode, n, psf, SBR, Polar, A, phi)
%   Note that Polar, A & phi are only needed in LDS modes
% -------------------------------------------------------------------------
% Coded by Zhengyi Zhan
% July 1, 2021
% -------------------------------------------------------------------------
narginchk(4,7);

if contains(mode, ["gaussian", "doughnut"], 'IgnoreCase', true)
    
    p0 = psf/sum(psf);
    p = (SBR*p0 + 1/length(psf)) / (SBR+1);
    out = -sum(n .* log(p));
    
else
    
    alpha = modeAlpha(mode);
    dotJ = A * cos(phi-alpha).^2 + (1-A)/2;
    numer = Polar .* psf .* dotJ;
    p0 = numer / sum(numer);
    p = (SBR*p0 + 1/length(psf)) / (SBR+1);
    out = -sum(n .* log(p));
    
end

end
