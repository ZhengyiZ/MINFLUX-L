%% PRE-PROCESS
addpath(genpath('Class'));
addpath(genpath('Function'));
addpath(genpath('Tool'));

clc; clear;

%% CALCULATION
% load datasets
load('Fig_4_a_to_c_Datasets.mat');
x0 = xStart:pixelSize:xStop;
y0 = yStart:pixelSize:yStop;

% calculate CRB
[psf,psfGx,psfGy] = numCal(mode,fwhm,xx,yy,L);
alpha = rad2deg(modeAlpha(mode));
for i = 1:length(x0)
    for j = 1:length(y0)
        [crbXY(i,j), crbA(i,j), crbPhi(i,j)] = ...
        covCalCRB(calCovMat(psf(i,j,:), psfGx(i,j,:), psfGy(i,j,:),...
        A0, alpha, phi0, SBR, nTotal));
    end
end

% RMSE and bias of Adam
[errXYAdam, biasXYAdam] = est2eb(xAdam, xx, yAdam, yy);
[errPhiAdam, biasPhiAdam] = est2eb4angle(rad2deg(phiAdam), phi0);
[errAAdam, biasAAdam] = est2eb(AAdam, A0);

% RMSE and bias of SGS
[errXYSGS, biasXYSGS] = est2eb(xSGS, xx, ySGS, yy);
[errPhiSGS, biasPhiSGS] = est2eb4angle(rad2deg(phiSGS), phi0);
[errASGS, biasASGS] = est2eb(ASGS, A0);

% RMSE and bias of mLMS
[errXYmLMS, biasXYmLMS] = est2eb(xmLMS, xx, ymLMS, yy);

%% DISPLAY
% Fig. 3(a)
figure,
imagesc(x0, y0, errXYAdam), axis image;
axis off; colorbar; title('RMSE_{xy} - Adam (nm)'); hold on;
set(gca,'clim',[1.5 5.5]);
% profile
dashline([-55;60], [0;0], 5, 2.5, 5, 2.5, 'color', 'w', 'LineWidth', 1);
% circle
circle_phi = 0:pi/40:2*pi;
circle_x = L/2*cos(circle_phi);
circle_y = L/2*sin(circle_phi);
plot(circle_x, circle_y, '--', 'Color', 'white', 'LineWidth', 1);
% PSF center position
SZPSF = 150;
LWPSF = 1;
scatter(0, 0, SZPSF, 'MarkerFaceColor', '#D95319',...
    'MarkerEdgeColor', '#FFFFFF', 'LineWidth', LWPSF);
scatter(0, L/2+0, SZPSF, 'MarkerFaceColor', '#A2142F',...
    'MarkerEdgeColor', '#FFFFFF', 'LineWidth', LWPSF);
scatter(0, 0-L/2, SZPSF, 'MarkerFaceColor', '#A2142F',...
    'MarkerEdgeColor', '#FFFFFF', 'LineWidth', LWPSF);
scatter(0+L/2*sind(60), 0+L/2*cosd(60), SZPSF, 'MarkerFaceColor', '#77AC30',...
    'MarkerEdgeColor', '#FFFFFF', 'LineWidth', LWPSF);
scatter(0-L/2*sind(60), 0-L/2*cosd(60), SZPSF, 'MarkerFaceColor', '#77AC30',...
    'MarkerEdgeColor', '#FFFFFF', 'LineWidth', LWPSF);
scatter(0+L/2*sind(120), 0+L/2*cosd(120), SZPSF, 'MarkerFaceColor', '#4DBEEE',...
    'MarkerEdgeColor', '#FFFFFF', 'LineWidth', LWPSF);
scatter(0-L/2*sind(120), 0-L/2*cosd(120), SZPSF, 'MarkerFaceColor', '#4DBEEE',...
    'MarkerEdgeColor', '#FFFFFF', 'LineWidth', LWPSF);
hold off;

% Fig. 3(b)
ind = 14;
figure,
plot(x0, crbXY(ind,:), 'LineWidth',2,'color','#EDB120'); hold on;
xlabel('Position (nm)'); ylabel('Precision_{xy} (nm)');
errorbar(x0, errXYAdam(ind,:), biasXYAdam(ind,:),'o','Color','#D95319',...
    'MarkerSize',5,'MarkerEdgeColor','#D95319','MarkerFaceColor','#D95319');
errorbar(x0, errXYSGS(ind,:), biasXYSGS(ind,:),'s','Color','#0072BD',...
    'MarkerSize',5,'MarkerEdgeColor','#0072BD','MarkerFaceColor','#0072BD');
errorbar(x0, errXYmLMS(ind,:), biasXYmLMS(ind,:),'d','Color','#7E2F8E',...
    'MarkerSize',5,'MarkerEdgeColor','#7E2F8E','MarkerFaceColor','#7E2F8E');
legend('CRB','Adam','GS','mLMS');
hold off;
set(gca,'ylim',[0 15]);
set(gca,'xlim',[-50 50]);

% Fig. 3(c)
figure,
subplot 211,
plot(x0, crbPhi(ind,:), 'LineWidth',2,'color','#EDB120'); hold on;
xlabel('Position (nm)'); ylabel('Precision_{\phi} (°)');
errorbar(x0, errPhiAdam(ind,:), biasPhiAdam(ind,:),'o','Color','#D95319',...
    'MarkerSize',5,'MarkerEdgeColor','#D95319','MarkerFaceColor','#D95319');
errorbar(x0, errPhiSGS(ind,:), biasPhiSGS(ind,:),'s','Color','#0072BD',...
    'MarkerSize',5,'MarkerEdgeColor','#0072BD','MarkerFaceColor','#0072BD');
legend('CRB','Adam','GS');
set(gca,'xlim',[-50 50]);
xticks([-50 -25 0 25 50]);
set(gca,'ylim',[3 8]);
subplot 212,
plot(x0, crbA(ind,:), 'LineWidth',2,'color','#EDB120'); hold on;
xlabel('Position (nm)'); ylabel('Precision_{A}');
errorbar(x0, errAAdam(ind,:), biasAAdam(ind,:),'o','Color','#D95319',...
    'MarkerSize',5,'MarkerEdgeColor','#D95319','MarkerFaceColor','#D95319');
errorbar(x0, errASGS(ind,:), biasASGS(ind,:),'s','Color','#0072BD',...
    'MarkerSize',5,'MarkerEdgeColor','#0072BD','MarkerFaceColor','#0072BD');
legend('CRB','Adam','GS');
set(gca,'xlim',[-50 50]);
xticks([-50 -25 0 25 50]);
set(gca,'ylim',[0.04 0.08]);
