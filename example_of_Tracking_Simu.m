%% PRE-PROCESS
addpath(genpath('Class'));
addpath(genpath('Function'));
addpath(genpath('Tool'));

clc; clear;

%% PARAMETERS
% simulation environment
fwhm            = 300;

% choose excitation beam pattern modes
mode            = 'LDS-3-9';

% beam center displacement
L               = 100;
SBR             = 3;
laserPower      = 150;
threshold       = 60;

% polarization non-uniformity, only valid for LDS dark spot
P               = [1 0.95 0.92];

% the number of tracking points
trackPoint      = 10000;

% lsc parameters
a               = 250;
b               = 250;
p               = 3;
q               = 2;
phi_lsc         = 0;
t               = linspace(0, 2*pi, trackPoint);
dipole.Offset.x = 0;
dipole.Offset.y = 0;

% polarization modulation depth
dipole.A0       = 0.5;

%% INITIALIZE
% create LMS object
lms = LMS(mode, fwhm, L, P);
lms.SBR = SBR;
lms.nVector = laserPower/10;
beta = lms.calBeta(1);

% calculate the Lissajous curve
[dipole.x, dipole.y, dipole.phi] = genLSC(a, b, p, q, phi_lsc, t);
dipole.x = dipole.x + dipole.Offset.x;
dipole.y = dipole.y + dipole.Offset.y;

%% LIVE TRAJECTORY
psf = numCal(mode, fwhm, 0, 0, L);
N = zeros(length(psf), trackPoint);
live.x = 0;
live.y = 0;
record.live.x = zeros(trackPoint, 1);
record.live.y = record.live.x;

textwaitbar(0, trackPoint/100, 'Tracking');

for i = 1:trackPoint

    % record current beam position before generate counts
    record.live.x(i) = live.x;
    record.live.y(i) = live.y;

    % calculate (x,y) point's psf & beam polarization
%     [p0, bg] = debyeGenPB(debye, dipole.x(i)-live.x, dipole.y(i)-live.y,...
%         0, L, SBR, dipole.A0, deg2rad(theta0), dipole.phi(i), P);
    [p0, bg] = numGenPB(mode, fwhm, dipole.x(i)-live.x, dipole.y(i)-live.y,...
        L, SBR, dipole.A0, deg2rad(dipole.phi(i)), P);
    n = phyGenCount1(p0, bg, laserPower);
    N(:,i) = n;

    % running lms
    lms.nVector = N(:,i);
    [xmLMS, ymLMS] = lms.calmLMS;

    % update live position only if the solution is reasonable
    if sqrt(xmLMS^2 + ymLMS^2) < threshold
        live.x = live.x + xmLMS;
        live.y = live.y + ymLMS;
    end

    if mod(i,100) == 0
        textwaitbar(i/100, trackPoint/100, 'Tracking');
    end

end

% N-hist
Nmean = mean(sum(N,1));
figure,
histogram(sum(N,1));
hold on;
plot([Nmean, Nmean], get(gca, 'YLim'), 'LineStyle','-', 'LineWidth', 2);
hold off;
fprintf('The average number of photons is %.2f.\n', Nmean);

% live trajectory
figure,
plot(record.live.x, record.live.y);
set(gca, 'xlim', [-300 300]);
set(gca, 'ylim', [-300 300]);
hold on;
patch(dipole.x, dipole.y, rad2deg(dipole.phi), 'EdgeColor', 'flat',...
    'FaceColor', 'none', 'LineWidth', 3);
colormap('hsv');
colorbar;
axis square;

%% POST-PROCESS
MLEobj = MLE(mode, fwhm, L, P);
MLEobj.SBR = SBR;
TWB = Timerwaitbar(trackPoint, 'Running Adam...');
for i = 1:trackPoint
    lms.nVector = N(:,i);
    [xmLMS, ymLMS] = lms.calmLMS;
    MLEobj.nVector = N(:,i);
    [xSolve, ySolve, ASolve, phiSolve] = MLEobj.Adam(0, xmLMS, ymLMS);
%     [xSolve, ySolve, ASolve, phiSolve] = MLEobj.GridSearch(0);
    if sqrt(xSolve^2 + ySolve^2) <= threshold
        record.post.x(i) = record.live.x(i) + xSolve;
        record.post.y(i) = record.live.y(i) + ySolve;
        record.post.A(i) = ASolve;
        record.post.phi(i) = phiSolve;
    else
        record.post.x(i) = NaN;
        record.post.y(i) = NaN;
        record.post.A(i) = NaN;
        record.post.phi(i) = NaN;
    end
    TWB.update();
end
    
TWB.delete();

%% DISPLAY
% errXY histogram & gamma fit
errXY = sqrt( ( (record.post.x-dipole.x).^2 +...
    (record.post.y-dipole.y).^2 ) );
errXY = errXY(~isnan(errXY));
figure, histogram(errXY);
[phat, pci] = gamfit(errXY);
fprintf('XY:\tμ = %.2f\tσ = %.2f\n', phat(1)*phat(2),...
    sqrt(phat(1))*phat(2));

if contains(mode, 'LDS')

    % errPhi histogram & norm fit
    errPhi = est2errV4angle(rad2deg(record.post.phi), rad2deg(dipole.phi));
    errPhi = errPhi(~isnan(errPhi));
    figure, histogram(errPhi);
    [muPhi, sigmaPhi] = normfit(errPhi);
    fprintf('Φ:\tμ = %.2f\tσ = %.2f\n', muPhi, sigmaPhi);

    % errA histogram & norm fit
    errA = record.post.A - dipole.A0;
    errA = errA(~isnan(errA));
    figure, histogram(errA);
    [muA, sigmaA] = normfit(errA);
    fprintf('A:\tμ = %.2f\tσ = %.2f\n', muA, sigmaA);

    % phi trajectory
    figure,
    patch(record.post.x, record.post.y, rad2deg(record.post.phi),...
        'EdgeColor', 'flat', 'FaceColor', 'none');
    colormap('hsv'); colorbar; axis square;
    set(gca, 'xlim', [-300 300]);
    set(gca, 'ylim', [-300 300]);
    
    % A trajectory
    figure,
    patch(record.post.x, record.post.y, record.post.A,...
        'EdgeColor', 'flat', 'FaceColor', 'none');
    colormap('default'); colorbar; axis square;
    set(gca, 'xlim', [-300 300]);
    set(gca, 'ylim', [-300 300]);

else

    % post trajectory
    figure,
    plot(record.post.x, record.post.y); axis square;
    set(gca, 'xlim', [-300 300]);
    set(gca, 'ylim', [-300 300]);

end
