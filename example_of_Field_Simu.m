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
mode            = 'doughnut-4';
% beam center displacement
L               = 100;
% total photon counts
nTotal          = 500;
% signal to background ratio
SBR             = 10;

% dipole position
x0              = -50:50;
y0              = x0;
pixelSize       = 5;

% polarization parameters, only valid for LDS
P               = [1 0.95 0.92];    % non-uniformity
A0              = 1;                % polarization modulation depth
phi0            = 75;               % azimuth angle, unit: degree.

laserPower      = 100;
timesEveryPoint = 10^1;

% Adam parameters
iterations      = 500;
lr              = 0.5;
lrPenalty       = 8;

%% INITIALIZE
xV = x0(1):pixelSize:x0(end);
yV = y0(1):pixelSize:y0(end);
[xx, yy] = meshgrid(xV, yV);

% generate p0 & bg
textwaitbar(0, length(yV), 'GEN P0 & BG');
for i = 1:length(yV)
    for j = 1:length(xV)
        [p0(i,j,:), bg(i,j)] = numGenPB(mode, fwhm, xx(i,j), yy(i,j), L,...
            SBR, A0, deg2rad(phi0), P);
    end
    textwaitbar(i, length(yV), 'GEN P0 & BG');
end

%% GENERATE PHOTON COUNTS
totalNum = numel(xx);
N = cell(length(yV),length(xV),timesEveryPoint);
bgTmp = zeros(length(yV),length(xV),timesEveryPoint);
SBRrecord = bgTmp;
Q = bgTmp;
for i = 1:length(yV)
    for j = 1:length(xV)
        for o = 1:timesEveryPoint
            [n, SBRest, count] = phyGenCount(p0(i,j,:), bg(i,j), nTotal, laserPower);
            N{i,j,o} = n;
            SBRrecord(i,j,o) = SBRest;
            Q(i,j,o) = count;
        end
        textwaitbar((i-1)*length(xV)+j, totalNum, 'GEN PHOTONS');
    end
end

%% LMS
% create the LMS object
LMSobj = LMS(mode, fwhm, L, P);
xLMS = zeros(length(yV), length(xV), timesEveryPoint);
yLMS = xLMS;
xmLMS = xLMS;
ymLMS = xLMS;

for i = 1:length(yV)
    for j = 1:length(xV)
        for o = 1:timesEveryPoint
            LMSobj.nVector = N{i,j,o};
            LMSobj.SBR = SBRrecord(i,j,o);
            [xLMS(i,j,o), yLMS(i,j,o)] = LMSobj.calLMS;
            if contains(mode, ["doughnut-4" "doughnut-7" "gaussian-4"...
                    "gaussian-7" "LDS-2-6" "LDS-3-9"], 'IgnoreCase', true)
                if ~exist('beta', 'var')
                    beta = LMSobj.calBeta(1);
                end
                [xmLMS(i,j,o), ymLMS(i,j,o)] = LMSobj.calmLMS;
            else
                xmLMS(i,j,o) = xLMS(i,j,o);
                ymLMS(i,j,o) = yLMS(i,j,o);
            end
        end
        textwaitbar((i-1)*length(xV)+j, totalNum, 'SOLVING LMS');
    end
end

%% MLE
% MLEobj = MLE(mode, fwhm, L, P);
% TWB = Timerwaitbar(totalNum, 'MLE');
% 
% for i = 1:length(yV)
%     for j = 1:length(xV)
%         for o = 1:timesEveryPoint
%             MLEobj.nVector = N{i,j,o};
%             MLEobj.SBR = SBRrecord(i,j,o);
%             [xAdam(i,j,o), yAdam(i,j,o), AAdam(i,j,o), phiAdam(i,j,o)] = ...
%                 MLEobj.Adam(0, 0, 0, iterations, lr, lrPenalty);
%             [xSGS(i,j,o), ySGS(i,j,o), ASGS(i,j,o), phiSGS(i,j,o)] = ...
%                 MLEobj.GridSearch(0);
%         end
%         TWB.update();
%     end
% end
% TWB.delete();

%% MLE USING PARALLEL POOL FOR ACCELERATION
TWB = Timerwaitbar(totalNum, 'MLE');

for i = 1:length(yV)
    for j = 1:length(xV)
        parfor o = 1:timesEveryPoint
            MLEobj = MLE(mode, fwhm, L, P);
            MLEobj.nVector = N{i,j,o};
            MLEobj.SBR = SBRrecord(i,j,o);
            [xAdam(i,j,o), yAdam(i,j,o), AAdam(i,j,o), phiAdam(i,j,o)] = ...
                MLEobj.Adam(0, 0, 0, iterations, lr, lrPenalty);
            [xSGS(i,j,o), ySGS(i,j,o), ASGS(i,j,o), phiSGS(i,j,o)] = ...
                MLEobj.GridSearch(0);
        end
        TWB.update();
    end
end
TWB.delete();

%% DISPLAY
crbXY = numCalCRB('xy', mode, fwhm, L, nTotal, SBR, xx, yy, A0, deg2rad(phi0));
[errXYAdam, biasXYAdam] = est2eb(xAdam, xx, yAdam, yy);
[errXYSGS, biasXYSGS] = est2eb(xSGS, xx, ySGS, yy);
[errXYmLMS, biasXYmLMS] = est2eb(xmLMS, xx, ymLMS, yy);
if contains(mode, 'LDS', 'IgnoreCase', true)
    crbA = numCalCRB('A', mode, fwhm, L, nTotal, SBR, xx, yy, A0, deg2rad(phi0));
    crbPhi = rad2deg(numCalCRB('phi', mode, fwhm, L, nTotal, SBR, xx, yy, A0, deg2rad(phi0)));
%     phiAdam = asin( sin(phiAdam) );
    [errPhiAdam, biasPhiAdam] = est2eb(rad2deg(phiAdam), phi0);
    [errAAdam, biasAAdam] = est2eb(AAdam, A0);
    [errPhiSGS, biasPhiSGS] = est2eb(rad2deg(phiSGS), phi0);
    [errASGS, biasASGS] = est2eb(ASGS, A0);
end

figure,
subplot 221, imagesc(xV, yV, crbXY), axis image;
colorbar; title('xy-CRB');
subplot 222, imagesc(xV, yV, errXYAdam), axis image;
colorbar; title('xy-Adam');
subplot 223, imagesc(xV, yV, errXYSGS), axis image;
colorbar; title('xy-SGS');
subplot 224, imagesc(xV, yV, errXYmLMS), axis image;
colorbar;
if contains(mode, ["doughnut-4" "doughnut-7" "gaussian-4"...
        "gaussian-7" "LDS-2-6" "LDS-3-9"], 'IgnoreCase', true)
    title('xy-mLMS');
else
    title('xy-LMS');
end

if contains(mode, 'LDS', 'IgnoreCase', true)
    figure,
    subplot 231, imagesc(xV, yV, crbA), axis image;
    colorbar; title('A-CRB');
    subplot 232, imagesc(xV, yV, errAAdam), axis image;
    colorbar; title('A-Adam');
    subplot 233, imagesc(xV, yV, errASGS), axis image;
    colorbar; title('A-SGS');
    subplot 234, imagesc(xV, yV, crbPhi), axis image;
    colorbar; title('\phi-CRB');
    subplot 235, imagesc(xV, yV, errPhiAdam), axis image;
    colorbar; title('\phi-Adam');
    subplot 236, imagesc(xV, yV, errPhiSGS), axis image;
    colorbar; title('\phi-SGS');
end
