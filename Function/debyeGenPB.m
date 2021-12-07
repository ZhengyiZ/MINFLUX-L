function [p0, bg] = debyeGenPB(debyeCal, x0, y0, z0, L, SBR, A0, theta0, phi0, P)
%DEBYEGENPB generates signal and background using debyeCal
% [p0, bg] = debyeGenPB(debyeCal, x0, y0, z0, L, SBR, A0, theta, phi, P)
% Note that A0, theta, phi & P are only necessary for 1D modes
% -------------------------------------------------------------------------
% Coded by Zhengyi Zhan
% July 1, 2021
% -------------------------------------------------------------------------
narginchk(6,10);

[exciPSF, orient] = debyeCal.calSP(x0, y0, z0, L);
if contains(debyeCal.mode, ["gaussian", "doughnut"], 'IgnoreCase', true)
    lambda = exciPSF;
else
    obj_orient = [sin(theta0) * cos(phi0)
        sin(theta0) * sin(phi0)
        cos(theta0)];
    lambda = exciPSF .* (A0 * dot(orient,...
        dotDimExpan(obj_orient, length(orient)), 2).^2 + (1-A0)/2);
    lambda = PExtend(P, debyeCal.mode) .* lambda;
end

p0 = lambda / sum(lambda);
bg = sum(p0) / SBR / length(p0);

end
