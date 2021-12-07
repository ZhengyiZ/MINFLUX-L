function [p0, bg] = numGenPB(mode, fwhm, x0, y0, L, SBR, A0, phi0, P)
%NUMGENPB generates signal and background using numCal
%   [p0, bg] = numGenPB(mode, fwhm, x0, y0, L, SBR, A0, phi0, P)
% Note that A0, phi0 & P are only necessary for LDS modes
% -------------------------------------------------------------------------
% Coded by Zhengyi Zhan
% July 1, 2021
% -------------------------------------------------------------------------
narginchk(6,9);

psf = numCal(mode, fwhm, x0, y0, L);
if contains(mode, ["gaussian", "doughnut"], 'IgnoreCase', true)
    lambda = psf;
else
    alpha = modeAlpha(mode);
    lambda = psf .* (A0 * cos(alpha-phi0).^2 + (1-A0)/2);
    if contains(mode, 'LDS-2')
        P = P(1:2);
    end
    lambda = PExtend(P, mode) .* lambda;
end

p0 = lambda / sum(lambda);
bg = sum(p0) / SBR / length(p0);

end
