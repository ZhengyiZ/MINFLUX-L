function fwhm = fminconFWHM(mode, exciWL, psf, x0, y0, alpha)
%FMINCONFWHM finds the FWHM of psf by mode
%   fwhm = findFWHM(mode, exciWL, psf, x, y, theta)
%   Note that alpha is only needed in LDS modes
% -------------------------------------------------------------------------
% Coded by Zhengyi Zhan
% July 1, 2021
% -------------------------------------------------------------------------
narginchk(5,6);

if length(x0) == numel(x0)
    [x,y] = meshgrid(x0, y0);
else
    x = x0;
    y = y0;
end

options = optimoptions(@fmincon, 'Display', 'off');
psfNorm = normCustom(psf);

if contains(mode, 'LDS', 'IgnoreCase', true)
    fun = @(fwhm) sum( abs( psfNorm - ...
        4*exp(1)*log(2)/fwhm^2 * (-x*sin(alpha)+y*cos(alpha)).^2 .*...
        exp(-4*log(2)/fwhm^2* (x.^2+y.^2) ) ) , 'all');
elseif contains(mode, 'doughnut', 'IgnoreCase', true)
    fun = @(fwhm) sum( abs( psfNorm - ...
        4*exp(1)*log(2)/fwhm^2 * (x.^2+y.^2) .* ...
        exp(-4*log(2)/fwhm^2* (x.^2+y.^2) ) ) , 'all');
elseif contains(mode, 'gaussian', 'IgnoreCase', true)
    fun = @(fwhm) sum( abs( psfNorm - ...
        exp(-4*log(2)/fwhm^2* (x.^2+y.^2) ) ) , 'all');
else
    error('Mode error');
end

fwhm = fmincon(fun,exciWL/3,0.5,exciWL,[],[],[],[],[],options);

end
