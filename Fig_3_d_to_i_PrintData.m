%% PRE-PROCESS
addpath(genpath('Class'));
addpath(genpath('Function'));
addpath(genpath('Tool'));

clc; clear;

%% CHOOSE DATASETS
flag = 'f';
switch flag
    case 'd'
        load('Fig_3_d_g_Datasets.mat');
        figTitle = 'y = 0';
    case 'e'
        load('Fig_3_e_h_Datasets.mat');
        figTitle = 'y = L/4';
    case 'f'
        load('Fig_3_f_i_Datasets.mat');
        figTitle = 'y = L/2';
end

%% CALCULATION
% initialization 
crbXY = zeros(length(L),length(nTotal));
errXYAdam = crbXY; biasXYAdam = crbXY;
errXYSGS = crbXY; biasXYSGS = crbXY;
crbA = zeros(size(nTotal));
crbPhi = zeros(size(nTotal));

% XY
for i = 1:length(L)
    crbXY(i,:) = numCalCRB('xy', mode, fwhm, L(i), nTotal, SBR,...
        x0(i), y0, A0, deg2rad(phi0));
    [errXYAdam(i,:), biasXYAdam(i,:)] = est2eb(xAdam(i,:,:), x0(i),...
        yAdam(i,:,:), y0);
    [errXYSGS(i,:), biasXYSGS(i,:)] = est2eb(xSGS(i,:,:), x0(i),...
        ySGS(i,:,:), y0);
end
% A
[errAAdam, biasAAdam] = est2eb(AAdam, A0);
[errASGS, biasASGS] = est2eb(ASGS, A0);
for i = 1:length(L)
    for j = 1:length(nTotal)
    crbA(i,j) = numCalCRB('A', mode, fwhm, L(i), nTotal(j), SBR,...
        x0(i), y0, A0, deg2rad(phi0));
    end
end

% Phi
[errPhiAdam, biasPhiAdam] = est2eb(rad2deg(phiAdam), phi0);
[errPhiSGS, biasPhiSGS] = est2eb(rad2deg(phiSGS), phi0);
for i = 1:length(L)
    for j = 1:length(nTotal)
        crbPhi(i,j) = rad2deg(numCalCRB('phi', mode, fwhm, L(i), nTotal(j),...
            SBR, x0(i), y0, A0, deg2rad(phi0)));
    end
end

%% DISPLAY
labels = {'L = 25 nm', 'L = 50 nm', 'L = 75 nm', 'L = 100 nm',...
    'L = 125 nm','L = 150 nm', 'Adam', 'SGS'};

% Fig. 3(d)/(e)/(f)
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
xlabel('N'); ylabel('error/nm');
legend([s1;s2;s3], labels, 'Location', 'best', 'NumColumns', 2);
title(figTitle);

figure,
subplot 211,
s21 = loglog(nTotal, crbPhi, 'color', '#A2142F');
hold on;
s22 = scatter(nTotal, errPhiAdam(1,:), 'MarkerEdgeColor', '#0072BD');
scatter(nTotal, errPhiAdam(2,:), 'MarkerEdgeColor', '#D95319');
scatter(nTotal, errPhiAdam(3,:), 'MarkerEdgeColor', '#EDB120');
scatter(nTotal, errPhiAdam(4,:), 'MarkerEdgeColor', '#7E2F8E');
scatter(nTotal, errPhiAdam(5,:), 'MarkerEdgeColor', '#77AC30');
scatter(nTotal, errPhiAdam(6,:), 'MarkerEdgeColor', '#4DBEEE');
s23 = scatter(nTotal, errPhiSGS(1,:), '*', 'MarkerEdgeColor', '#0072BD');
scatter(nTotal, errPhiSGS(2,:), '*', 'MarkerEdgeColor', '#D95319');	
scatter(nTotal, errPhiSGS(3,:), '*', 'MarkerEdgeColor', '#EDB120');
scatter(nTotal, errPhiSGS(4,:), '*', 'MarkerEdgeColor', '#7E2F8E');	
scatter(nTotal, errPhiSGS(5,:), '*', 'MarkerEdgeColor', '#77AC30');
scatter(nTotal, errPhiSGS(6,:), '*', 'MarkerEdgeColor', '#4DBEEE');
hold off;
set(gca, 'xlim', [3 1300]);
legend([s21(1), s22, s23], "CRB", "Adam", "SGS", 'Location', 'best');

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
suptitle(figTitle);
