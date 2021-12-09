%% PRE-PROCESS
addpath(genpath('Class'));
addpath(genpath('Function'));
addpath(genpath('Tool'));

clc;
clear;
close all;

%% PARAMETERS SETTING
% simulation environment
fwhm            = 300;

% choose excitation modes from 
%             'xxxxxxxx-3' | 'xxxxxxxx-4' | 'xxxxxxxx-6' | 'xxxxxxxx-7'
%             'LDS-2-4' | 'LDS-2-5' | 'LDS-2-6'
%             'LDS-3-6' | 'LDS-3-7' | 'LDS-3-9'
%             'xxxxxxxx' could be 'Gaussian' or 'Doughnut'
mode            = 'LDS-3-9';
% beam center displacement
L               = 100;
% total photon counts
N               = 100;
% signal to background ratio
SBR             = 10;

x0              = -300:300;
y0              = x0;
pixelSize       = 4;

% polarization parameters
A0              = 0.5;        % the value range of A should be [0,1].
phi0            = 30;         % azimuth angle, unit: degree.

%% CALCULATE CRB FOR XY IN XY PLANE
% display PSF
xV = x0(1):pixelSize:x0(end);
yV = y0(1):pixelSize:y0(end);
psf = numCal(mode, fwhm, xV, yV, L);
dispPSF(psf, mode);

% calculate CRB
crbXY = numCalCRB('xy', mode, fwhm, L, N, SBR, x0, y0, A0, deg2rad(phi0));
figure, imagesc(xV, yV, crbXY); axis image; colorbar;
title('Precision_{xy} (nm)');
crbXYmin = min(crbXY, [], 'all');
set(gca, 'clim', [crbXYmin crbXYmin+50]);
% profile
figure, plot( crbXY(:, round(size(crbXY,2)/2)) );
xlabel('x Position (nm)'); ylabel('Precision_{xy} (nm)');

%% CALCULATE CRB FOR A/PHI IN XY PLANE
if contains(mode, 'LDS', 'IgnoreCase', true)
    crbA = numCalCRB('A', mode, fwhm, L, N, SBR, x0, y0, A0, deg2rad(phi0));
    crbPhi = rad2deg(numCalCRB('phi', mode, fwhm, L, N, SBR, x0, y0, A0, deg2rad(phi0)));
    figure,
    subplot 121, imagesc(xV, yV, crbA); axis image; colorbar;
    title('Precision_A');
    subplot 122, imagesc(xV, yV, crbPhi); axis image; colorbar;
    title('Precision_\phi (°）');
else
    disp('Only LDS modes can calculate polarization parameters.');
end

%% CALCULATE CRB FOR XY ALONG N
x0 = 0;
y0 = 0;
N = unique(round(logspace(0.4,4.1,100)));
crbXY_N = numCalCRB('xy', mode, fwhm, L, N, SBR, x0, y0, A0, deg2rad(phi0));
figure,
loglog(N, crbXY_N); axis tight;
xlabel('Total number of photons'); ylabel('Precision_{xy} (nm)');

%% CALCULTING CRB FOR A/PHI ALONG N
if contains(mode, 'LDS', 'IgnoreCase', true)
    crbA_N = numCalCRB('A', mode, fwhm, L, N, SBR, x0, y0, A0, deg2rad(phi0));
    crbPhi_N = rad2deg(numCalCRB('phi', mode, fwhm, L, N, SBR, x0, y0, A0, deg2rad(phi0)));
    figure,
    subplot 121, loglog(N, crbA_N); axis tight;
    xlabel('Total number of photons'); ylabel('Precision_{A}');
    subplot 122, loglog(N, crbPhi_N); axis tight;
    xlabel('Total number of photons'); ylabel('Precision_{\phi} (°)');
end

%% CALCULTING CRB FOR A/PHI IN A-PHI PLANE
if contains(mode, 'LDS', 'IgnoreCase', true)
    N = 100;
    A = 0:0.01:1;
    phi = deg2rad(-90:90);
    [AA, pp] = meshgrid(A, phi);
    for i = 1:length(phi)
        for j = 1:length(A)
            crbA_ap(i,j) = numCalCRB('A', mode, fwhm, L, N, SBR, x0, y0,...
                AA(i,j), pp(i,j));
            crbPhi_ap(i,j) = rad2deg(numCalCRB('phi', mode, fwhm, L, N,...
                SBR, x0, y0, AA(i,j), pp(i,j)));
        end
    end
    figure,
    subplot 121, imagesc(A, phi, crbA_ap); axis square; colorbar;
    title('Precision_{A}'); set(gca, 'clim', [0 0.2]);
    subplot 122, imagesc(A, phi, crbPhi_ap); axis square; colorbar;
    title('Precision_{\phi} (°)'); set(gca, 'clim', [0 45])
end
