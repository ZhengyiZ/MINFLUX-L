%% PRE-PROCESS
addpath(genpath('Class'));
addpath(genpath('Function'));
addpath(genpath('Tool'));

clc; clear;
close all;

colorA = '#A1D6E4';
colorB = '#FFa67a';
colorC = '#CCFFCC';
colorAS = '#2f8ca5';
colorBS = '#EF5000';
colorCS = '#00a500';

%% LOAD DATASETS
load('Fig_5_a_15.mat');
NA = N;
recordA = record;

load('Fig_5_b_31.mat');
NB = N;
recordB = record;

load('Fig_5_c_48.mat');
NC = N;
recordC = record;

%% CALCULATION
% number of photons
TNA = sum(NA, 1);
TNB = sum(NB, 1);
TNC = sum(NC, 1);
muNA = poissfit(TNA);
muNB = poissfit(TNB);
muNC = poissfit(TNC);

% XY error
errXYA = sqrt( ( (recordA.post.x-dipole.x).^2 +...
    (recordA.post.y-dipole.y).^2 ) );
errXYB = sqrt( ( (recordB.post.x-dipole.x).^2 +...
    (recordB.post.y-dipole.y).^2 ) );
errXYC = sqrt( ( (recordC.post.x-dipole.x).^2 +...
    (recordC.post.y-dipole.y).^2 ) );
errXYA = errXYA(~isnan(errXYA));
errXYB = errXYB(~isnan(errXYB));
errXYC = errXYC(~isnan(errXYC));
phatA = gamfit(errXYA); muXYA = prod(phatA); sigmaXYA = sqrt(phatA(1))*phatA(2);
phatB = gamfit(errXYB); muXYB = prod(phatB); sigmaXYB = sqrt(phatB(1))*phatB(2);
phatC = gamfit(errXYC); muXYC = prod(phatC); sigmaXYC = sqrt(phatC(1))*phatC(2);

% Phi error
errPhiA = est2errV4angle(rad2deg(recordA.post.phi), rad2deg(dipole.phi));
errPhiB = est2errV4angle(rad2deg(recordB.post.phi), rad2deg(dipole.phi));
errPhiC = est2errV4angle(rad2deg(recordC.post.phi), rad2deg(dipole.phi));
errPhiA = errPhiA(~isnan(errPhiA));
errPhiB = errPhiB(~isnan(errPhiB));
errPhiC = errPhiC(~isnan(errPhiC));
phiDistA = fitdist(errPhiA', 'tLocationScale');
phiDistB = fitdist(errPhiB', 'tLocationScale');
phiDistC = fitdist(errPhiC', 'tLocationScale');

% A error
errAA = recordA.post.A - dipole.A0;
errAB = recordB.post.A - dipole.A0;
errAC = recordC.post.A - dipole.A0;
errAA = errAA(~isnan(errAA));
errAB = errAB(~isnan(errAB));
errAC = errAC(~isnan(errAC));
[muAA, sigmaAA] = normfit(errAA);
[muAB, sigmaAB] = normfit(errAB);
[muAC, sigmaAC] = normfit(errAC);

%% DISPLAY TRAJECTORY
figure,
subplot 231,
plot(recordA.live.x, recordA.live.y);
set(gca, 'xlim', [-300 300]); set(gca, 'ylim', [-300 300]);
hold on;
patch(dipole.x, dipole.y, rad2deg(dipole.phi), 'EdgeColor', 'flat',...
    'FaceColor', 'none', 'LineWidth', 2);
colormap('hsv'); axis square;
colorbar('Ticks',[-90,-45,0,45,90]);

subplot 232,
plot(recordB.live.x, recordB.live.y);
set(gca, 'xlim', [-300 300]); set(gca, 'ylim', [-300 300]);
hold on;
patch(dipole.x, dipole.y, rad2deg(dipole.phi), 'EdgeColor', 'flat',...
    'FaceColor', 'none', 'LineWidth', 2);
colormap('hsv'); axis square;
colorbar('Ticks',[-90,-45,0,45,90]);

subplot 233,
plot(recordC.live.x, recordC.live.y);
set(gca, 'xlim', [-300 300]); set(gca, 'ylim', [-300 300]);
hold on;
patch(dipole.x, dipole.y, rad2deg(dipole.phi), 'EdgeColor', 'flat',...
    'FaceColor', 'none', 'LineWidth', 2);
colormap('hsv'); axis square;
colorbar('Ticks',[-90,-45,0,45,90]);

subplot 234,
patch(recordA.post.x, recordA.post.y, rad2deg(recordA.post.phi),...
    'EdgeColor', 'flat', 'FaceColor', 'none');
set(gca, 'xlim', [-300 300]); set(gca, 'ylim', [-300 300]);
colormap('hsv'); axis square;
colorbar('Ticks',[-90,-45,0,45,90]);

subplot 235,
patch(recordB.post.x, recordB.post.y, rad2deg(recordB.post.phi),...
    'EdgeColor', 'flat', 'FaceColor', 'none');
set(gca, 'xlim', [-300 300]); set(gca, 'ylim', [-300 300]);
colormap('hsv'); axis square;
colorbar('Ticks',[-90,-45,0,45,90]);

subplot 236,
patch(recordC.post.x, recordC.post.y, rad2deg(recordC.post.phi),...
    'EdgeColor', 'flat', 'FaceColor', 'none');
set(gca, 'xlim', [-300 300]); set(gca, 'ylim', [-300 300]);
colormap('hsv'); axis square;
colorbar('Ticks',[-90,-45,0,45,90]);

figure,
subplot 131,
patch(recordA.post.x, recordA.post.y, recordA.post.A, ...
    'EdgeColor', 'flat', 'FaceColor', 'none');
set(gca, 'xlim', [-300 300]); set(gca, 'ylim', [-300 300]);
colormap('default'); axis square;
colorbar('Ticks',[0 0.25 0.5 0.75 1]);

subplot 132,
patch(recordB.post.x, recordB.post.y, recordB.post.A, ...
    'EdgeColor', 'flat', 'FaceColor', 'none');
set(gca, 'xlim', [-300 300]); set(gca, 'ylim', [-300 300]);
colormap('default'); axis square;
colorbar('Ticks',[0 0.25 0.5 0.75 1]);

subplot 133,
patch(recordC.post.x, recordC.post.y, recordC.post.A, ...
    'EdgeColor', 'flat', 'FaceColor', 'none');
set(gca, 'xlim', [-300 300]); set(gca, 'ylim', [-300 300]);
colormap('default'); axis square;
colorbar('Ticks',[0 0.25 0.5 0.75 1]);

%% DISPLAY HISTOGRAMS
% number of photons
figure,
histogram( TNA ,'Normalization', 'pdf', 'BinWidth', 1, ...
    'FaceColor', colorA);
hold on;
histogram( TNB ,'Normalization', 'pdf', 'BinWidth', 1, ...
    'FaceColor', colorB);
histogram( TNC ,'Normalization', 'pdf', 'BinWidth', 1, ...
    'FaceColor', colorC);

set(gca, 'XLim', [0 70]);
xN = 0:70;
yNA = poisspdf(xN, muNA);
yNB = poisspdf(xN, muNB);
yNC = poisspdf(xN, muNC);

plot(xN, yNA, 'LineWidth', 2, 'Color', colorAS);
plot(xN, yNB, 'LineWidth', 2, 'Color', colorBS);
plot(xN, yNC, 'LineWidth', 2, 'Color', colorCS);
hold off;

legend(['<N> = ' num2str(muNA,'%.2f')], ['<N> = ' num2str(muNB,'%.2f')], ...
    ['<N> = ' num2str(muNC,'%.2f')]);
xlabel('Total number of photons');
ylabel('Probabilty Density');

figure,
% xy error
subplot(4,3,[1 4]),
histogram( errXYA ,'Normalization', 'pdf', 'BinWidth', 1.5, ...
    'FaceColor', colorA, 'FaceAlpha', 0.6);
hold on;
set(gca, 'XLim', [0 50]);
set(gca, 'YLim', [0 0.1]);
xA = linspace(0, 50, 1000);
yXYA = gampdf(xA, phatA(1), phatA(2));
plot(xA, yXYA, 'LineWidth', 2, 'Color', colorAS);
hold off;
xlabel('Error_{xy} (nm)');
ylabel('Probabilty Density');

subplot(4,3,[2 5]),
histogram( errXYB ,'Normalization', 'pdf', 'BinWidth', 1.5, ...
    'FaceColor', colorB, 'FaceAlpha', 0.6);
hold on;
set(gca, 'XLim', [0 50]);
set(gca, 'YLim', [0 0.1]);
yXYB = gampdf(xA, phatB(1), phatB(2));
plot(xA, yXYB, 'LineWidth', 2, 'Color', colorBS);
hold off;
xlabel('Error_{xy} (nm)');
ylabel('Probabilty Density');

subplot(4,3,[3 6]),
histogram( errXYC ,'Normalization', 'pdf', 'BinWidth', 1.5, ...
    'FaceColor', colorC, 'FaceAlpha', 0.6);
hold on;
set(gca, 'XLim', [0 50]);
set(gca, 'YLim', [0 0.1]);
yXYC = gampdf(xA, phatC(1), phatC(2));
plot(xA, yXYC, 'LineWidth', 2, 'Color', colorCS);
hold off;
xlabel('Error_{xy} (nm)');
ylabel('Probabilty Density');

% phi error
subplot(4,3,7),
histogram( errPhiA ,'Normalization', 'pdf', 'BinWidth', 5, ...
    'FaceColor', colorA, 'FaceAlpha', 0.6);
hold on;
set(gca, 'XLim', [-90 90]);
set(gca, 'YLim', [0 0.03]);
xA = linspace(-90, 90, 1000);
yPhiA = phiDistA.pdf(xA);
plot(xA, yPhiA, 'LineWidth', 2, 'Color', colorAS);
hold off;
xlabel('Error_{\Phi} (°)');
ylabel('Probability Density');

subplot(4,3,8),
histogram( errPhiB ,'Normalization', 'pdf', 'BinWidth', 5, ...
    'FaceColor', colorB, 'FaceAlpha', 0.6);
hold on;
set(gca, 'XLim', [-90 90]);
set(gca, 'YLim', [0 0.03]);
yPhiB = phiDistB.pdf(xA);
plot(xA, yPhiB, 'LineWidth', 2, 'Color', colorBS);
hold off;
xlabel('Error_{\Phi} (°)');
ylabel('Probability Density');

subplot(4,3,9),
histogram( errPhiC ,'Normalization', 'pdf', 'BinWidth', 5, ...
    'FaceColor', colorC, 'FaceAlpha', 0.6);
hold on;
set(gca, 'XLim', [-90 90]);
set(gca, 'YLim', [0 0.03]);
yPhiC = phiDistC.pdf(xA);
plot(xA, yPhiC, 'LineWidth', 2, 'Color', colorCS);
hold off;
xlabel('Error_{\Phi} (°)');
ylabel('Probability Density');

% A error
subplot(4,3,10),
histogram( errAA ,'Normalization', 'pdf', 'BinWidth', 0.025, ...
    'FaceColor', colorA, 'FaceAlpha', 0.6);
hold on;
set(gca, 'XLim', [-0.5 0.5]);
set(gca, 'YLim', [0 7]);
xA = linspace(-0.5, 0.5, 1000);
yAA = normpdf(xA, muAA, sigmaAA);
plot(xA, yAA, 'LineWidth', 2, 'Color', colorAS);
hold off;
xlabel('Error_{A}');
ylabel('Probability Density');

subplot(4,3,11),
histogram( errAB ,'Normalization', 'pdf', 'BinWidth', 0.025, ...
    'FaceColor', colorB, 'FaceAlpha', 0.6);
hold on;
set(gca, 'XLim', [-0.5 0.5]);
set(gca, 'YLim', [0 7]);
yAB = normpdf(xA, muAB, sigmaAB);
plot(xA, yAB, 'LineWidth', 2, 'Color', colorBS);
hold off;
xlabel('Error_{A}');
ylabel('Probability Density');

subplot(4,3,12),
histogram( errAC ,'Normalization', 'pdf', 'BinWidth', 0.025, ...
    'FaceColor', colorC, 'FaceAlpha', 0.6);
hold on;
set(gca, 'XLim', [-0.5 0.5]);
set(gca, 'YLim', [0 7]);
yAC = normpdf(xA, muAC, sigmaAC);
plot(xA, yAC, 'LineWidth', 2, 'Color', colorCS);
hold off;
xlabel('Error_{A}');
ylabel('Probability Density');
