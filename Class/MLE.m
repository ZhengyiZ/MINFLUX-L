%% MLE USING ADAM AND GRID SEARCH ALGORITHMS
%  Version    : v2.0
%  Author     : Zhengyi, Zhan
%  Release    : 9th Sep. 2021
%  Dependency : numCal.m
%
%  USAGE
%  Creation of MLE
%  - MLEobj = MLE(mode, fwhm, L, P);
%    INPUT PARAMETERS
%       mode: 'xxxxxxxx-3' | 'xxxxxxxx-4' | 'xxxxxxxx-6' | 'xxxxxxxx-7'
%             'LDS-2-4' | 'LDS-2-5' | 'LDS-2-6'
%             'LDS-3-6' | 'LDS-3-7' | 'LDS-3-9'
%             'xxxxxxxx' could be 'Gaussian' or 'Doughnut'
%       fwhm: full width at half height, unit: nm
%       L: scaning pattern displacement distance, unit: nm
%       P: polarization non-uniformity,
%          not necessary for Gaussian & Doughnut modes
%
%  Set SBR - signal to background ratio
%  - MLEobj.SBR = SBR;
%  Set nVector
%  - MLEobj.nVector = nVector;
%  Set P
%  - MLEobj.P = P;
%
%  Method for Adam
%  - [x, y, A, phi] = MLEobj.Adam(visFlag, iterations, lr, lrPenalty,...
%                            beta1, beta2, eps);
%    INPUT PARAMETERS
%       visualFlag: 1 for enable figure, 0 for disable, default: 0
%       iterations: iterations of Adam, recommended value depends on
%                   different L, such as 500 for L = 100 nm
%                   default: 500
%       lr: learing rate of Adam, default: 0.5
%       lrPenalty: learing rate penalty for A & phi, default: 10
%       Note that the following parameters are not recommended to be changed
%       beta1: default: 0.9
%       beta2: default: 0.999
%       eps: default: 1e-8
%    OUTPUT
%       x: estimated value of x
%       y: estimated value of y
%       A: estimated value of a, invalid in gaussian and doughnut modes
%       phi: estimated value of dipole orientation in xy-plane
%  - [x, y, A, phi] = MLEobj.Adam(visFlag, xS, yS, iterations,...
%                            lr, lrPenalty, beta1, beta2, eps);
%    INPUT & OUTPUT PARAMETERS are as above
%       xS: start point of x
%       yS: start point of x
%
%  Method for Grid Search
%  - [x, y, A, phi] = GridSearch(visFlag);
%    INPUT & OUTPUT PARAMETERS are as above
%
%  EXAMPLE
%  % Adam parameters
%  iterations      = 500;
%  lr              = 0.5;
%  lrPenalty       = 10;
%  visFlag         = 1;
%
%  % Creation of an MLEobj object
%  MLEobj = MLEobj(mode, fwhm, L, P);
%
%  % Set the nVector & SBR to be estimated before running algorithms
%  MLEobj.nVector = n;
%  MLEobj.SBR = SBR;
%
%  % Run Adam algorithm
%  [x, y, A, phi] = MLEobj.Adam(visFlag, 0, 0, iterations, lr, lrPenalty);
%
%  % Run Grid Search algorithm
%  [x, y, A, phi] = MLEobj.GridSearch(visFlag);
%
%  Public Properties
%  - MLEobj.SBR
%  - MLEobj.nVector
%  - MLEobj.P

classdef MLE < handle
    properties (Access = private)
        mode
        L
        fwhm
    end
    
    properties (Access = public)
        P
        SBR
        nVector
    end
    
    methods
        
        function obj = MLE(mode, fwhm, L, P)
            
            narginchk(3,4);
            obj.mode = mode;
            obj.fwhm = fwhm;
            obj.L = L;
            
            if contains(mode, 'LDS', 'IgnoreCase', true)
                obj.P = P;
            end
            
        end
        
        function [x, y, A, phi] = Adam(obj, visFlag, xS, yS, iterations,...
            lr, lrPenalty, beta1, beta2, eps)
            
            if isempty(obj.nVector)
                error('Please set nVector before running Adam.');
            elseif isempty(obj.SBR)
                error('Please set SBR before running Adam.');
            elseif isempty(obj.P)
                if contains(obj.mode, 'LDS', 'IgnoreCase', true)
                    error('Please set P before running Adam.');
                end
            end
            
            if nargin < 5
                iterations = 500;
            end
            
            if nargin < 7
                lr = 1;
                lrPenalty = 100;
            end
            
            if nargin < 10
                beta1 = 0.9;
                beta2 = 0.999;
                eps = 1e-8;
            end
            
            A = 0.5; phi = 0; x = xS; y = yS;
            m = 0; v = 0;
            if visFlag
                eval = zeros(iterations+1, 1);
            end
            
            if contains(obj.mode, ["gaussian", "doughnut"], 'IgnoreCase', true)
                
                Polar = 1;
                
                if visFlag
                    History = zeros(iterations,2);
                    equal_lr = zeros(iterations,2);
                end
                
                for t = 1:iterations
                    
                    % calculate intensity & gradient simutaneously
                    [psf, psfGx, psfGy] = numCal(obj.mode, obj.fwhm, x, y, obj.L);
                    psf(psf == 0) = 1e-18;
                    
                    % calculate current evaluation function
                    if visFlag
                        eval(t) = evalFun(obj.mode, obj.nVector, psf, obj.SBR);
                    end
                    
                    psfSum = sum(psf);
                    p = psf ./ psfSum ;
                    dpdx = ( psfSum.*psfGx - psf.*sum(psfGx ) ) ./ psfSum.^2;
                    dpdy = ( psfSum.*psfGy - psf.*sum(psfGy) ) ./ psfSum.^2;
                    grad = -[sum( obj.nVector .* (dpdx./(1/length(psf)/obj.SBR+p)) ) ...
                        sum( obj.nVector .* (dpdy./(1/length(psf)/obj.SBR+p)) )];
                    
                    % Adam routine
                    m = beta1 * m + (1-beta1) * grad;
%                     v = beta2 * v + (1-beta2) * grad.^2;
                    v = max(beta2 * v + (1-beta2) * grad.^2, v);
                    m_tip = m / (1-beta1^t);
                    v_tip = v / (1-beta2^t);
                    delta = - lr .* m_tip ./ (sqrt(v_tip)+eps);
                    if visFlag
                        equal_lr(t,:) = - delta ./ grad;
                    end
                    
                    % Adam update
                    x = x + delta(1);
                    y = y + delta(2);
                    
                    % record history
                    if visFlag
                        History(t,:) = [x;y];
                    end
                    
                end
                
            else
                
                alpha = modeAlpha(obj.mode);
                Polar = PExtend(obj.P, obj.mode);
                lr = [ones(1,2)*lr ones(1,1)*lr/lrPenalty ...
                    ones(1,1)*lr/lrPenalty*2*pi];
                if visFlag
                    equal_lr = zeros(iterations,4);
                    History = zeros(iterations,4);
                end
                
                for t = 1:iterations
                    
                    % calculate intensity & gradient simutaneously
                    [psf, psfGx, psfGy] = numCal(obj.mode, obj.fwhm, x, y, obj.L);
                    
                    % calculate current evaluation function
                    if visFlag
                        eval(t) = evalFun(obj.mode, obj.nVector, psf,...
                            obj.SBR, Polar, A, phi);
                    end
                    
                    polarItem = A * cos(phi-alpha).^2 + (1-A)/2;
                    u = Polar .* psf .* polarItem;
                    su = sum(u);
                    p = u ./ su;
                    dudx = Polar .* polarItem .* psfGx;
                    dudy = Polar .* polarItem .* psfGy;
                    dudA = Polar .* psf .* cos(2*(alpha-phi)) /2;
                    dudPhi = A * Polar .* psf .* sin(2*(alpha-phi));
                    dpdx = ( su.*dudx - u.*sum(dudx) ) ./ su.^2;
                    dpdy = ( su.*dudy - u.*sum(dudy) ) ./ su.^2;
                    dpdA = ( su.*dudA - u.*sum(dudA) ) ./ su.^2;
                    dpdAlpha = ( su.*dudPhi - u.*sum(dudPhi) ) ./ su.^2;
                    
                    % calculate gradient vector
                    grad = -[sum(obj.nVector .* (dpdx./(1/length(psf)/obj.SBR+p)) ) ...
                        sum(obj.nVector .* (dpdy./(1/length(psf)/obj.SBR+p)) ) ...
                        sum(obj.nVector .* (dpdA./(1/length(psf)/obj.SBR+p)) )...
                        sum(obj.nVector .* (dpdAlpha./(1/length(psf)/obj.SBR+p)) )];
                    
                    % Adam routine
                    m = beta1 * m + (1-beta1) * grad;
%                     v = beta2 * v + (1-beta2) * grad.^2;
                    v = max(beta2 * v + (1-beta2) * grad.^2, v);
                    m_tip = m / (1-beta1^t);
                    v_tip = v / (1-beta2^t);
                    delta = - lr .* m_tip ./ (sqrt(v_tip)+eps);
                    % Nadam
%                     delta = - lr .* (beta1*m_tip+(1-beta1)/(1-beta1^t)*grad) ./ (sqrt(v_tip)+eps);
                    if visFlag
                        equal_lr(t,:) = - delta ./ grad;
                    end
                    
                    % Adam update
                    x = x + delta(1);
                    y = y + delta(2);
                    A = A + delta(3);
                    A = min(max(A,1e-5),1);         % limit the range of At
                    phi = asin(sin(phi + delta(4)));
                    
                    % record history
                    if visFlag
                        History(t,:) = [x;y;A;rad2deg(phi)];
                    end
                    
                end
            end
            
            if visFlag
                
                psfTmp = numCal(obj.mode, obj.fwhm, x, y, obj.L);
                eval(t+1) = evalFun(obj.mode, obj.nVector, psfTmp, obj.SBR, Polar, A, phi);
                
                if contains(obj.mode, ["gaussian", "doughnut"], 'IgnoreCase', true)
                    figure, sgtitle(['lr=' num2str(lr(1)) ', T=' num2str(iterations)...
                        ', \beta_1=' num2str(beta1) ', \beta_2=' num2str(beta2)...
                        ', \epsilon=' num2str(eps)]);
                    subplot 221, plot(eval), legend('eval');
                    subplot 223,  plot(equal_lr(:,1));
                    hold on, plot(equal_lr(:,2));
                    hold off; legend('lr_{equal}X','lr_{equal}Y');
                    subplot 222, plot(History(:,1)), ylabel('x');
                    subplot 224, plot(History(:,2)), ylabel('y');
                else
                    figure, sgtitle(['lr=' num2str(lr(1)) ', lrPenalty=' num2str(lrPenalty)...
                        ', T=' num2str(iterations) ', \beta_1=' num2str(beta1) ', \beta_2=' num2str(beta2)...
                        ', \epsilon=' num2str(eps)]);
                    subplot(4,4,[1 2 5 6]), plot(eval), legend('eval');
                    subplot(4,4,[9 10 13 14]), plot(equal_lr(:,1));
                    hold on, plot(equal_lr(:,2));
                    hold on, plot(equal_lr(:,3));
                    hold on, plot(equal_lr(:,4));
                    hold off; legend('lr_{equal}X','lr_{equal}Y','lr_{equal}A',...
                        'lr_{equal}\phi');
                    subplot(4,4,[3 4]), plot(History(:,1)), ylabel('x');
                    subplot(4,4,[7 8]), plot(History(:,2)), ylabel('y');
                    subplot(4,4,[11 12]), plot(History(:,3)), ylabel('A');
                    subplot(4,4,[15 16]), plot(History(:,4)), ylabel('\phi');
                end
     
            end
            
        end
        
        function [x, y, A, phi] = GridSearch(obj, visFlag)
            
            if isempty(obj.nVector)
                error('Please set nVector before running Grid Search.');
            elseif isempty(obj.SBR)
                error('Please set SBR before running Grid Search.');
            elseif isempty(obj.P)
                if contains(obj.mode, 'LDS', 'IgnoreCase', true)
                    error('Please set P before running Grid Search.');
                end
            end
            
            if nargin < 2
                visFlag = 0;
            end
            
            if contains(obj.mode, ["gaussian", "doughnut"], 'IgnoreCase', true)
                
                A = 0;
                phi = 0;
                [x0,y0] = meshgrid(-200:5:200,-200:5:200);
                [x,y,res] = searchWorker(obj, x0, y0);
                if visFlag
                    figure,
                    subplot 221, imagesc(x0(1,:),y0(:,1),res); axis image;
                    hold on;
                    plot(x,y,'r.','markersize',10);
                    hold off;
                end
                
                [x0,y0] = meshgrid(x-5:1:x+5,y-5:1:y+5);
                [x,y,res] = searchWorker(obj, x0, y0);
                if visFlag
                    subplot 222, imagesc(x0(1,:),y0(:,1),res); axis image;
                    hold on;
                    plot(x,y,'r.','markersize',10);
                    hold off;
                end
                
                [x0,y0] = meshgrid(x-0.5:0.1:x+0.5,y-0.5:0.1:y+0.5);
                [x,y,res] = searchWorker(obj, x0, y0);
                if visFlag
                    subplot 223, imagesc(x0(1,:),y0(:,1),res); axis image;
                    hold on;
                    plot(x,y,'r.','markersize',10);
                    hold off;
                end
                
                [x0,y0] = meshgrid(x-0.1:0.01:x+0.1,y-0.1:0.01:y+0.1);
                [x,y,res] = searchWorker(obj, x0, y0);
                if visFlag
                    subplot 224, imagesc(x0(1,:),y0(:,1),res); axis image;
                    hold on;
                    plot(x,y,'r.','markersize',10);
                    hold off;
                end
                
            elseif contains(obj.mode, 'LDS', 'IgnoreCase', true)
                
                [x0,y0] = meshgrid(-100:5:100,-100:5:100);
                A0 = 0:0.1:1;
                phi0 = deg2rad(-90:10:90);
                [x, y, A, phi] = searchWorkerLDS(obj, x0, y0, A0, phi0);
                
                [x0,y0] = meshgrid(x-5:1:x+5,y-5:1:y+5);
                A0 = max(A-0.25,0):0.05:min(A+0.25,1);
                phi = rad2deg(phi);
                phi0 = deg2rad(max(phi-45,-90):2:min(phi+45,90));
                [x, y, A, phi] = searchWorkerLDS(obj, x0, y0, A0, phi0);
                
                [x0,y0] = meshgrid(x-1:0.1:x+1,y-1:0.1:y+1);
                A0 = max(A-0.1,0):0.01:min(A+0.1,1);
                phi = rad2deg(phi);
                phi0 = deg2rad(max(phi-20,-90):0.5:min(phi+20,90));
                [x, y, A, phi] = searchWorkerLDS(obj, x0, y0, A0, phi0);
     
            end

        end
        
    end
    
end

function [x, y, res] = searchWorker(obj, x0, y0)

psf = numCal(obj.mode,obj.fwhm,x0,y0,obj.L);
if obj.SBR == Inf
    p = psf ./ sum(psf,3);
else
    p = (obj.SBR*psf./sum(psf,3) + 1/size(psf,3)) / (obj.SBR+1);
end
res = -sum(dotDimExpan(obj.nVector) .* log(p),3);
x = x0(res == min(res,[],'all'));
y = y0(res == min(res,[],'all'));
if length(x) > 1
    x = x(1);
end
if length(y) > 1
    y = y(1);
end

end

function [x, y, A, phi] = searchWorkerLDS(obj, x0, y0, A0, phi0)

alpha = modeAlpha(obj.mode);
Polar = PExtend(obj.P, obj.mode);
res = zeros([size(x0) length(A0) length(phi0)]);
psf = numCal(obj.mode,obj.fwhm,x0,y0,obj.L);
for i = 1:length(A0)
    for j = 1:length(phi0)
        psfDot = psf .* dotDimExpan(Polar .* (A0(i)*cos(phi0(j)-alpha).^2+(1-A0(i))/2));
        if obj.SBR == Inf
            p = psfDot ./ sum(psfDot,3);
        else
            p = (obj.SBR*psfDot./sum(psfDot,3) + 1/size(psfDot,3)) / (obj.SBR+1);
        end
        res(:,:,i,j) = -sum(dotDimExpan(obj.nVector) .* log(p),3);
    end
end
ind = find(res == min(res,[],'all'));
[sx,sy,sA,sa] = ind2sub(size(res),ind);
x = x0(sx,sy);
y = y0(sx,sy);
A = A0(sA);
phi = phi0(sa);
if length(x) > 1
    x = x(1);
end
if length(y) > 1
    y = y(1);
end
if length(A) > 1
    A = A(1);
end
if length(phi) > 1
    phi = phi(1);
end

end
                