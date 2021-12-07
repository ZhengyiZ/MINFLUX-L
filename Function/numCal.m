function [psf, psfGx, psfGy] = numCal(mode, fwhm, x, y, L)
%NUMCAL PSF calculation using the approximate expression
%   [psf, psfGx, psfGy] = numCal(mode, fwhm, x, y, L)
%   calculates the psf and the gradient of psf
%
%   INPUT PARAMETERS
%      mode: 'xxxxxxxx-3' | 'xxxxxxxx-4' | 'xxxxxxxx-6' | 'xxxxxxxx-7'
%            'LDS-2-4' | 'LDS-2-5' | 'LDS-2-6'
%            'LDS-3-6' | 'LDS-3-7' | 'LDS-3-9'     
%            'xxxxxxxx' could be 'Gaussian' or 'Doughnut'
%      L: scaning pattern displacement distance, unit: nm
% -------------------------------------------------------------------------
% Coded by Zhengyi Zhan
% Sept 9, 2021
% -------------------------------------------------------------------------

if isvector(x) && isvector(y)
    [x,y] = meshgrid(x,y);
elseif isvector(x) || isvector(y)
    error('x y should be the same type, vectors or matrixes.');
end

if contains(mode, ["gaussian", "doughnut"], 'IgnoreCase', true)
    if contains(mode, ["3", "4"])
        psf = zeros([size(x) 3]);
        psfGx = zeros([size(x) 3]);
        psfGy = zeros([size(x) 3]);
        for i = 1:3
            if contains(mode, 'doughnut', 'IgnoreCase', true)
                psf(:,:,i) = doughnutPSF(fwhm, x+L/2*cosd(120*i),...
                    y+L/2*sind(120*i));
                psfGx(:,:,i) = doughnutGx(fwhm, x+L/2*cosd(120*i),...
                    y+L/2*sind(120*i));
                psfGy(:,:,i) = doughnutGy(fwhm, x+L/2*cosd(120*i),...
                    y+L/2*sind(120*i));
            else
                psf(:,:,i) = gaussianPSF(fwhm, x+L/2*cosd(120*i),...
                    y+L/2*sind(120*i));
                psfGx(:,:,i) = gaussianGx(fwhm, x+L/2*cosd(120*i),...
                    y+L/2*sind(120*i));
                psfGy(:,:,i) = gaussianGy(fwhm, x+L/2*cosd(120*i),...
                    y+L/2*sind(120*i));
            end
        end
        if contains(mode, '4')
            if contains(mode, 'doughnut', 'IgnoreCase', true)
                psf(:,:,i+1) = doughnutPSF(fwhm, x, y);
                psfGx(:,:,i+1) = doughnutGx(fwhm, x, y);
                psfGy(:,:,i+1) = doughnutGy(fwhm, x, y);
            else
                psf(:,:,i+1) = gaussianPSF(fwhm, x, y);
                psfGx(:,:,i+1) = gaussianGx(fwhm, x, y);
                psfGy(:,:,i+1) = gaussianGy(fwhm, x, y);
            end
        end
    elseif contains(mode, ["6", "7"])
        psf = zeros([size(x) 6]);
        psfGx = zeros([size(x) 6]);
        psfGy = zeros([size(x) 6]);
        for i = 1:6
            if contains(mode, 'doughnut', 'IgnoreCase', true)
                psf(:,:,i) = doughnutPSF(fwhm, x+L/2*cosd(60*i),...
                    y+L/2*sind(60*i));
                psfGx(:,:,i) = doughnutGx(fwhm, x+L/2*cosd(60*i),...
                    y+L/2*sind(60*i));
                psfGy(:,:,i) = doughnutGy(fwhm, x+L/2*cosd(60*i),...
                    y+L/2*sind(60*i));
            else
                psf(:,:,i) = gaussianPSF(fwhm, x+L/2*cosd(60*i),...
                    y+L/2*sind(60*i));
                psfGx(:,:,i) = gaussianGx(fwhm, x+L/2*cosd(60*i),...
                    y+L/2*sind(60*i));
                psfGy(:,:,i) = gaussianGy(fwhm, x+L/2*cosd(60*i),...
                    y+L/2*sind(60*i));
            end
        end
        if contains(mode, '7')
            if contains(mode, 'doughnut', 'IgnoreCase', true)
                psf(:,:,i+1) = doughnutPSF(fwhm, x, y);
                psfGx(:,:,i+1) = doughnutGx(fwhm, x, y);
                psfGy(:,:,i+1) = doughnutGy(fwhm, x, y);
            else
                psf(:,:,i+1) = gaussianPSF(fwhm, x, y);
                psfGx(:,:,i+1) = gaussianGx(fwhm, x, y);
                psfGy(:,:,i+1) = gaussianGy(fwhm, x, y);
            end
        end
    else
        error('mode error.');
    end
    
elseif contains(mode, 'LDS-2', 'IgnoreCase', true)
    
    phi = [0;90];
    
    if contains(mode, '4')
        psf = zeros([size(x) 4]);
        psfGx = zeros([size(x) 4]);
        psfGy = zeros([size(x) 4]);
        for i = 1:length(phi)
            psf(:,:,2*i-1) = LDSPSF(fwhm, x-L/2*sind(phi(i)),...
                y+L/2*cosd(phi(i)), phi(i));
            psf(:,:,2*i) = LDSPSF(fwhm, x+L/2*sind(phi(i)),...
                y-L/2*cosd(phi(i)), phi(i));
            psfGx(:,:,2*i-1) = LDSGx(fwhm, x-L/2*sind(phi(i)),...
                y+L/2*cosd(phi(i)), phi(i));
            psfGx(:,:,2*i) = LDSGx(fwhm, x+L/2*sind(phi(i)),...
                y-L/2*cosd(phi(i)), phi(i));
            psfGy(:,:,2*i-1) = LDSGy(fwhm, x-L/2*sind(phi(i)),...
                y+L/2*cosd(phi(i)), phi(i));
            psfGy(:,:,2*i) = LDSGy(fwhm, x+L/2*sind(phi(i)),...
                y-L/2*cosd(phi(i)), phi(i));
        end
    elseif contains(mode, '5')
        psf = zeros([size(x) 5]);
        psfGx = zeros([size(x) 5]);
        psfGy = zeros([size(x) 5]);
        for i = 1:length(phi)
            psf(:,:,3*i-2) = LDSPSF(fwhm, x-L/2*sind(phi(i)),...
                y+L/2*cosd(phi(i)), phi(i));
            psf(:,:,3*i-1) = LDSPSF(fwhm, x+L/2*sind(phi(i)),...
                y-L/2*cosd(phi(i)), phi(i));
            psfGx(:,:,3*i-2) = LDSGx(fwhm, x-L/2*sind(phi(i)),...
                y+L/2*cosd(phi(i)), phi(i));
            psfGx(:,:,3*i-1) = LDSGx(fwhm, x+L/2*sind(phi(i)),...
                y-L/2*cosd(phi(i)), phi(i));
            psfGy(:,:,3*i-2) = LDSGy(fwhm, x-L/2*sind(phi(i)),...
                y+L/2*cosd(phi(i)), phi(i));
            psfGy(:,:,3*i-1) = LDSGy(fwhm, x+L/2*sind(phi(i)),...
                y-L/2*cosd(phi(i)), phi(i));
        end
        psf(:,:,3) = LDSPSF(fwhm, x, y, phi(1));
        psfGx(:,:,3) = LDSGx(fwhm, x, y, phi(1));
        psfGy(:,:,3) = LDSGy(fwhm, x, y, phi(1));
    elseif contains(mode, '6')
        psf = zeros([size(x) 6]);
        psfGx = zeros([size(x) 6]);
        psfGy = zeros([size(x) 6]);
        for i = 1:length(phi)
            psf(:,:,3*i-2) = LDSPSF(fwhm, x-L/2*sind(phi(i)),...
                y+L/2*cosd(phi(i)), phi(i));
            psf(:,:,3*i-1) = LDSPSF(fwhm, x+L/2*sind(phi(i)),...
                y-L/2*cosd(phi(i)), phi(i));
            psfGx(:,:,3*i-2) = LDSGx(fwhm, x-L/2*sind(phi(i)),...
                y+L/2*cosd(phi(i)), phi(i));
            psfGx(:,:,3*i-1) = LDSGx(fwhm, x+L/2*sind(phi(i)),...
                y-L/2*cosd(phi(i)), phi(i));
            psfGy(:,:,3*i-2) = LDSGy(fwhm, x-L/2*sind(phi(i)),...
                y+L/2*cosd(phi(i)), phi(i));
            psfGy(:,:,3*i-1) = LDSGy(fwhm, x+L/2*sind(phi(i)),...
                y-L/2*cosd(phi(i)), phi(i));
            psf(:,:,3*i) = LDSPSF(fwhm, x, y, phi(i));
            psfGx(:,:,3*i) = LDSGx(fwhm, x, y, phi(i));
            psfGy(:,:,3*i) = LDSGy(fwhm, x, y, phi(i));
        end
    else
        error('mode error.');
    end
    
elseif contains(mode, 'LDS-3', 'IgnoreCase', true)
    
    phi = [0;60;-60];
    
    if contains(mode, '6')
        psf = zeros([size(x) 6]);
        psfGx = zeros([size(x) 6]);
        psfGy = zeros([size(x) 6]);
        for i = 1:length(phi)
            psf(:,:,2*i-1) = LDSPSF(fwhm, x-L/2*sind(phi(i)),...
                y+L/2*cosd(phi(i)), phi(i));
            psf(:,:,2*i) = LDSPSF(fwhm, x+L/2*sind(phi(i)),...
                y-L/2*cosd(phi(i)), phi(i));
            psfGx(:,:,2*i-1) = LDSGx(fwhm, x-L/2*sind(phi(i)),...
                y+L/2*cosd(phi(i)), phi(i));
            psfGx(:,:,2*i) = LDSGx(fwhm, x+L/2*sind(phi(i)),...
                y-L/2*cosd(phi(i)), phi(i));
            psfGy(:,:,2*i-1) = LDSGy(fwhm, x-L/2*sind(phi(i)),...
                y+L/2*cosd(phi(i)), phi(i));
            psfGy(:,:,2*i) = LDSGy(fwhm, x+L/2*sind(phi(i)),...
                y-L/2*cosd(phi(i)), phi(i));
        end
    elseif contains(mode, '7')
        psfTmp = zeros([size(x) 8]);
        psfGxTmp = zeros([size(x) 8]);
        psfGyTmp = zeros([size(x) 8]);
        psf = zeros([size(x) 7]);
        psfGx = zeros([size(x) 7]);
        psfGy = zeros([size(x) 7]);
        for i = 1:length(phi)
            psfTmp(:,:,3*i-2) = LDSPSF(fwhm, x-L/2*sind(phi(i)),...
                y+L/2*cosd(phi(i)), phi(i));
            psfTmp(:,:,3*i-1) = LDSPSF(fwhm, x+L/2*sind(phi(i)),...
                y-L/2*cosd(phi(i)), phi(i));
            psfGxTmp(:,:,3*i-2) = LDSGx(fwhm, x-L/2*sind(phi(i)),...
                y+L/2*cosd(phi(i)), phi(i));
            psfGxTmp(:,:,3*i-1) = LDSGx(fwhm, x+L/2*sind(phi(i)),...
                y-L/2*cosd(phi(i)), phi(i));
            psfGyTmp(:,:,3*i-2) = LDSGy(fwhm, x-L/2*sind(phi(i)),...
                y+L/2*cosd(phi(i)), phi(i));
            psfGyTmp(:,:,3*i-1) = LDSGy(fwhm, x+L/2*sind(phi(i)),...
                y-L/2*cosd(phi(i)), phi(i));
        end
        psf(:,:,[1:2,4:5,6:7]) = psfTmp(:,:,[1:2,4:5,7:8]);
        psfGx(:,:,[1:2,4:5,6:7]) = psfGxTmp(:,:,[1:2,4:5,7:8]);
        psfGy(:,:,[1:2,4:5,6:7]) = psfGyTmp(:,:,[1:2,4:5,7:8]);
        psf(:,:,3) = LDSPSF(fwhm, x, y, phi(1));
        psfGx(:,:,3) = LDSGx(fwhm, x, y, phi(1));
        psfGy(:,:,3) = LDSGy(fwhm, x, y, phi(1));
    elseif contains(mode, '9')
        psf = zeros([size(x) 9]);
        psfGx = zeros([size(x) 9]);
        psfGy = zeros([size(x) 9]);
        for i = 1:length(phi)
            psf(:,:,3*i-2) = LDSPSF(fwhm, x-L/2*sind(phi(i)),...
                y+L/2*cosd(phi(i)), phi(i));
            psf(:,:,3*i-1) = LDSPSF(fwhm, x+L/2*sind(phi(i)),...
                y-L/2*cosd(phi(i)), phi(i));
            psfGx(:,:,3*i-2) = LDSGx(fwhm, x-L/2*sind(phi(i)),...
                y+L/2*cosd(phi(i)), phi(i));
            psfGx(:,:,3*i-1) = LDSGx(fwhm, x+L/2*sind(phi(i)),...
                y-L/2*cosd(phi(i)), phi(i));
            psfGy(:,:,3*i-2) = LDSGy(fwhm, x-L/2*sind(phi(i)),...
                y+L/2*cosd(phi(i)), phi(i));
            psfGy(:,:,3*i-1) = LDSGy(fwhm, x+L/2*sind(phi(i)),...
                y-L/2*cosd(phi(i)), phi(i));
            psf(:,:,3*i) = LDSPSF(fwhm, x, y, phi(i));
            psfGx(:,:,3*i) = LDSGx(fwhm, x, y, phi(i));
            psfGy(:,:,3*i) = LDSGy(fwhm, x, y, phi(i));
        end
    else
        error('mode error.');
    end
else
    error('mode error.');
end

if ndims(psf) == 3 && numel(psf) == length(psf)
    psf = squeeze(psf);
    psfGx = squeeze(psfGx);
    psfGy = squeeze(psfGy);
end

end

function res = gaussianPSF(fwhm, x, y)
res = exp(-4*log(2)/fwhm^2* (x.^2+y.^2) );
end

function res = gaussianGx(fwhm, x, y)
res = -8*log(2)*x/fwhm^2.*exp(-4*log(2)/fwhm^2*(x.^2+y.^2));
end

function res = gaussianGy(fwhm, x, y)
res = -8*log(2)*y/fwhm^2.*exp(-4*log(2)/fwhm^2*(x.^2+y.^2));
end

function res = doughnutPSF(fwhm, x, y)
res = 4*exp(1)*log(2)/fwhm^2 * (x.^2+y.^2) .*...
    exp(-4*log(2)/fwhm^2* (x.^2+y.^2) );
end

function res = doughnutGx(fwhm, x, y)
res = 8*exp(1)*log(2)*x/fwhm^4.*(fwhm^2-4*log(2)*(x.^2+y.^2)).*...
    exp(-4*log(2)/fwhm^2*(x.^2+y.^2));
end

function res = doughnutGy(fwhm, x, y)
res = 8*exp(1)*log(2)*y/fwhm^4.*(fwhm^2-4*log(2)*(x.^2+y.^2)).*...
    exp(-4*log(2)/fwhm^2*(x.^2+y.^2));
end

function res = LDSPSF(fwhm, x, y, phi)
res = 4*exp(1)*log(2)/fwhm^2 * (-x*sind(phi)+y*cosd(phi)).^2 .*...
    exp(-4*log(2)/fwhm^2* (x.^2+y.^2) );
end

function res = LDSGx(fwhm, x, y, phi)
res = 8*exp(1)*log(2)/fwhm^4 * (x*sind(phi)-y*cosd(phi)) .* ...
    (4*log(2)*x.*y*cosd(phi)+sind(phi)*(fwhm^2-4*log(2)*x.^2)) .*...
    exp(-4*log(2)/fwhm^2* (x.^2+y.^2) );
end

function res = LDSGy(fwhm, x, y, phi)
res = -8*exp(1)*log(2)/fwhm^4 * (x*sind(phi)-y*cosd(phi)) .* ...
    (4*log(2)*x.*y*sind(phi)+cosd(phi)*(fwhm^2-4*log(2)*y.^2)) .*...
    exp(-4*log(2)/fwhm^2* (x.^2+y.^2) );
end
