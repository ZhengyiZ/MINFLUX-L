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
nTotal          = 13;
laserPower      = 60;
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
lms.nVector = nTotal;
beta = lms.calBeta(1);

% calculate the Lissajous curve
[dipole.x, dipole.y, dipole.phi] = genLSC(a, b, p, q, phi_lsc, t);
dipole.x = dipole.x + dipole.Offset.x;
dipole.y = dipole.y + dipole.Offset.y;

%% LIVE TRAJECTORY
psf = numCal(mode, fwhm, 0, 0, L);
N = zeros(length(psf), trackPoint);
Q = zeros(1, trackPoint);
live.x = 0;
live.y = 0;
record.live.x = zeros(trackPoint, 1);
record.live.y = record.live.x;
TWB = Timerwaitbar(trackPoint, 'Processing...');

for i = 1:trackPoint
        
        % record current beam position before generate counts
        record.live.x(i) = live.x;
        record.live.y(i) = live.y;
        
        dipole.orient = cos(dipole.phi(i));
        
        % calculate (x,y) point's psf & beam polarization
        % [psf, orient] = debye.calSP(x(i)-xCurr,y(i)-yCurr,0,L);
        psf = numCal(mode, fwhm, dipole.x(i)-live.x, dipole.y(i)-live.y, L);
        
        if contains(mode, ["gaussian", "doughnut"], 'IgnoreCase', true)
            lambda = psf;
        else
            % lambda = PExtend(P, mode) .* psf .* ...
            % (A * dot(orient,dotDimExpan(obj_orient,length(orient)),2).^2 + (1-A)/2);
            lambda = PExtend(P, mode) .* psf .* ...
                (dipole.A0 * cos(dipole.phi(i)-modeAlpha(mode)).^2 + (1-dipole.A0)/2);
        end
        p0 = lambda / sum(lambda) * (laserPower/10);
        bg = sum(p0) / (SBR+1) / length(p0);
        
        pData = zeros(nTotal*length(p0), length(p0));
        bgData = pData;
        for k = 1:length(p0)
            pData(:,k) = poissrnd(p0(k), nTotal*length(p0), 1);
            bgData(:,k) = poissrnd(bg, nTotal*length(p0), 1);
        end
        pb = pData + bgData;
        n = zeros(length(p0),1);
        count = 1;
        while sum(n) < nTotal
            if count > 1 && sum(n)+sum(pb(count,:)) > 1.5*nTotal
                break;
            end
            for k = 1:length(p0)
                n(k) = n(k) + pb(count, k);
            end
            count = count + 1;
            if count > length(pData)
                error('not enough.');
            end
        end
        Q(:,i) = count;
        N(:,i) = n;
        
        % running lms
        lms.nVector = N(:,i);
        [xmLMS, ymLMS] = lms.calmLMS;
        
        % update live position only if the solution is reasonable
        if sqrt(xmLMS^2 + ymLMS^2) < threshold
            live.x = live.x + xmLMS;
            live.y = live.y + ymLMS;
        end
    
    TWB.update();
    
end

TWB.delete();

figure,
plot(record.live.x,record.live.y);
set(gca, 'xlim', [-300 300]);
set(gca, 'ylim', [-300 300]);
hold on;
patch(dipole.x,dipole.y,dipole.phi,'EdgeColor','flat','FaceColor','none','LineWidth',3);
colormap('hsv');
colorbar;
axis square;

disp(['The average number of photons is ' num2str(mean(sum(N,1))) '.']);
disp(['The average number of EBP is ' num2str(mean(Q)) '.']);

%% POST-PROCESS
MLEobj = MLE(mode, fwhm, L, P);
MLEobj.SBR = SBR;
TWB = Timerwaitbar(trackPoint, 'Running Adam...');
for i = 1:trackPoint
    lms.nVector = N(:,i);
    [xmLMS, ymLMS] = lms.calmLMS;
    MLEobj.nVector = N(:,i);
    [xAdam, yAdam, AAdam, phiAdam] = MLEobj.Adam(0, xmLMS, ymLMS);
    if sqrt(xAdam^2 + yAdam^2) <= threshold
        record.post.x(i) = record.live.x(i) + xAdam;
        record.post.y(i) = record.live.y(i) + yAdam;
        record.post.A(i) = AAdam;
        record.post.phi(i) = phiAdam;
    else
        if i ~= 1
            record.post.x(i) = record.post.x(i-1);
            record.post.y(i) = record.post.y(i-1);
            record.post.A(i) = AAdam;
            record.post.phi(i) = phiAdam;
        else
            record.post.x(i) = 0;
            record.post.y(i) = 0;
            record.post.A(i) = AAdam;
            record.post.phi(i) = phiAdam;
        end
    end
    TWB.update();
end
    
TWB.delete();

if contains(mode, 'LDS')
    figure,
    patch(record.post.x, record.post.y, record.post.phi,...
        'EdgeColor','flat','FaceColor','none');
    colormap('hsv'); colorbar; axis square;
    set(gca, 'xlim', [-300 300]);
    set(gca, 'ylim', [-300 300]);
    
    figure,
    patch(record.post.x, record.post.y, record.post.A,...
        'EdgeColor','flat','FaceColor','none');
    colormap('default'); colorbar; axis square;
    set(gca, 'xlim', [-300 300]);
    set(gca, 'ylim', [-300 300]);
else
    figure,
    plot(record.post.x, record.post.y); axis square;
    set(gca, 'xlim', [-300 300]);
    set(gca, 'ylim', [-300 300]);
end
