%% PRE-PROCESS
addpath(genpath('Class'));
addpath(genpath('Function'));
addpath(genpath('Tool'));

clc;
clear;

%% PARAMETERS SETTING
% simulation environment
fwhm            = 300;

% choose excitation beam pattern modes
mode            = 'LDS-3-9';

L               = [25 50 75 100 125 150];
SBR             = 10;
nTotal          = unique(round(logspace(0.4,3.1,20)));

% dipole position
x0              = L/4;
y0              = 0;

% polarization parameters, only valid for LDS
P               = [1 0.95 0.92];    % non-uniformity
A0              = 1;                % polarization modulation depth
phi0            = 75;               % azimuth angle, unit: degree.

timesEveryPoint = 10^3;

% Adam parameters
iterations      = 1000;
lr              = 1;
lrPenalty       = 10;

%% INITIALIZE
N = cell(length(L), length(nTotal), timesEveryPoint);
totalNum = length(L) * length(nTotal);

for i = 1:length(L)
    
    % generate p0 & bg
    [p0(:,i), bg(i)] = numGenPB(mode, fwhm, x0(i), y0, L(i), SBR, A0, deg2rad(phi0), P);
    p(:,i) = (p0(:,i) + bg(i)) / sum(p0(:,i) + bg(i));

    % Generate photon counts
    pDist = makedist('Multinomial', 'probabilities', p(:,i));
    
    for j = 1:length(nTotal)
        for o = 1:timesEveryPoint
            n = distGenCount(pDist, nTotal(j));
            N{i,j,o} = n;
        end
    end
    
    textwaitbar(i, length(L), 'GEN PHOTONS');
    
end

%% MLE
xAdam = zeros(length(L), length(nTotal), timesEveryPoint);
yAdam = xAdam; AAdam = xAdam; phiAdam = xAdam;
xSGS = xAdam; ySGS = xAdam; ASGS = xAdam; phiSGS = xAdam;
TWB = Timerwaitbar(totalNum);
for i = 1:length(L)
    for j = 1:length(nTotal)
        parfor o = 1:timesEveryPoint
            MLEobj = MLE(mode, fwhm, L(i), P);
            MLEobj.nVector = N{i,j,o};
            MLEobj.SBR = SBR;
            [xAdam(i,j,o), yAdam(i,j,o), AAdam(i,j,o), phiAdam(i,j,o)] =...
                MLEobj.Adam(0, 0, 0, iterations, lr, lrPenalty);
            [xSGS(i,j,o), ySGS(i,j,o), ASGS(i,j,o), phiSGS(i,j,o)] =...
                MLEobj.GridSearch(0);
        end
        TWB.update();
    end
end
TWB.delete();

%% DISPLAY
crbXY = zeros([length(L) length(nTotal)]);
errXYAdam = crbXY; biasXYAdam = crbXY;
errXYSGS = crbXY; biasXYSGS = crbXY;
for i = 1:length(L)
    crbXY(i,:) = numCalCRB('xy', mode, fwhm, L(i), nTotal, SBR, x0(i), y0, A0, deg2rad(phi0));
    [errXYAdam(i,:), biasXYAdam(i,:)] = est2eb(xAdam(i,:,:), x0(i), yAdam(i,:,:), y0);
    [errXYSGS(i,:), biasXYSGS(i,:)] = est2eb(xSGS(i,:,:), x0(i), ySGS(i,:,:), y0);
end

figure,
s1 = loglog(nTotal, crbXY);
hold on;
s2 = scatter(nTotal, errXYAdam(1,:), 'MarkerEdgeColor', '#0072BD');
scatter(nTotal, errXYAdam(2,:), 'MarkerEdgeColor', '#D95319');	
scatter(nTotal, errXYAdam(3,:), 'MarkerEdgeColor', '#EDB120');
scatter(nTotal, errXYAdam(4,:), 'MarkerEdgeColor', '#7E2F8E');	
scatter(nTotal, errXYAdam(5,:), 'MarkerEdgeColor', '#77AC30');
scatter(nTotal, errXYAdam(6,:), 'MarkerEdgeColor', '#4DBEEE');
s3 = scatter(nTotal, errXYSGS(1,:), '*', 'MarkerEdgeColor', '#0072BD');
scatter(nTotal, errXYSGS(2,:), '*', 'MarkerEdgeColor', '#D95319');	
scatter(nTotal, errXYSGS(3,:), '*', 'MarkerEdgeColor', '#EDB120');
scatter(nTotal, errXYSGS(4,:), '*', 'MarkerEdgeColor', '#7E2F8E');	
scatter(nTotal, errXYSGS(5,:), '*', 'MarkerEdgeColor', '#77AC30');
scatter(nTotal, errXYSGS(6,:), '*', 'MarkerEdgeColor', '#4DBEEE');
hold off;
set(gca, 'xlim', [3 1300]);
set(gca, 'ylim', [0.1 200]);
xlabel('N'); ylabel('RMSE_{xy} (nm)');
legend([s1;s2;s3], "L = 25 nm", "L = 50 nm", "L = 75 nm", "L = 100 nm",...
    "L = 125 nm", "L = 150 nm", "Adam", "SGS", 'Location', 'best');

if contains(mode, 'LDS', 'IgnoreCase', true)
    crbA = crbXY;
    crbPhi = crbA;
    [errAAdam, biasAAdam] = est2eb(AAdam, A0);
    [errASGS, biasASGS] = est2eb(ASGS, A0);
    [errPhiAdam, biasPhiAdam] = est2eb4angle(rad2deg(phiAdam), phi0);
    [errPhiSGS, biasPhiSGS] = est2eb4angle(rad2deg(phiSGS), phi0);
    for i = 1:length(L)
        for j = 1:length(nTotal)
            crbA(i,j) = numCalCRB('A', mode, fwhm, L(i), nTotal(j), SBR,...
                x0(i), y0, A0, deg2rad(phi0));
            crbPhi(i,j) = rad2deg(numCalCRB('phi', mode, fwhm, L(i), nTotal(j),...
                SBR, x0(i), y0, A0, deg2rad(phi0)));
        end
    end
    
    figure,
    subplot 211,
    s1 = loglog(nTotal, crbPhi, 'color', '#A2142F');
    hold on;
    s2 = scatter(nTotal, errPhiAdam(1,:), 'MarkerEdgeColor', '#0072BD');
    scatter(nTotal, errPhiAdam(2,:), 'MarkerEdgeColor', '#D95319');
    scatter(nTotal, errPhiAdam(3,:), 'MarkerEdgeColor', '#EDB120');
    scatter(nTotal, errPhiAdam(4,:), 'MarkerEdgeColor', '#7E2F8E');
    scatter(nTotal, errPhiAdam(5,:), 'MarkerEdgeColor', '#77AC30');
    scatter(nTotal, errPhiAdam(6,:), 'MarkerEdgeColor', '#4DBEEE');
    s3 = scatter(nTotal, errPhiSGS(1,:), '*', 'MarkerEdgeColor', '#0072BD');
    scatter(nTotal, errPhiSGS(2,:), '*', 'MarkerEdgeColor', '#D95319');
    scatter(nTotal, errPhiSGS(3,:), '*', 'MarkerEdgeColor', '#EDB120');
    scatter(nTotal, errPhiSGS(4,:), '*', 'MarkerEdgeColor', '#7E2F8E');
    scatter(nTotal, errPhiSGS(5,:), '*', 'MarkerEdgeColor', '#77AC30');
    scatter(nTotal, errPhiSGS(6,:), '*', 'MarkerEdgeColor', '#4DBEEE');
    hold off;
    set(gca, 'xlim', [3 1300]);
    xlabel('N'); ylabel('RMSE_{\phi} (°)');
    legend([s1(1);s2;s3], "CRB", "Adam", "SGS", 'Location', 'best');
    
    subplot 212, loglog(nTotal, crbA, 'color', '#A2142F');
    hold on;
    scatter(nTotal, errAAdam(1,:), 'MarkerEdgeColor', '#0072BD');
    scatter(nTotal, errAAdam(2,:), 'MarkerEdgeColor', '#D95319');
    scatter(nTotal, errAAdam(3,:), 'MarkerEdgeColor', '#EDB120');
    scatter(nTotal, errAAdam(4,:), 'MarkerEdgeColor', '#7E2F8E');
    scatter(nTotal, errAAdam(5,:), 'MarkerEdgeColor', '#77AC30');
    scatter(nTotal, errAAdam(6,:), 'MarkerEdgeColor', '#4DBEEE');
    scatter(nTotal, errASGS(1,:), '*', 'MarkerEdgeColor', '#0072BD');
    scatter(nTotal, errASGS(2,:), '*', 'MarkerEdgeColor', '#D95319');
    scatter(nTotal, errASGS(3,:), '*', 'MarkerEdgeColor', '#EDB120');
    scatter(nTotal, errASGS(4,:), '*', 'MarkerEdgeColor', '#7E2F8E');
    scatter(nTotal, errASGS(5,:), '*', 'MarkerEdgeColor', '#77AC30');
    scatter(nTotal, errASGS(6,:), '*', 'MarkerEdgeColor', '#4DBEEE');
    hold off;
    set(gca, 'xlim', [3 1300]);
    xlabel('N'); ylabel('RMSE_{A}');
    legend([s1(1);s2;s3], "CRB", "Adam", "SGS", 'Location', 'best');
end
