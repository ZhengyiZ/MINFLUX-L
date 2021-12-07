%% LMS LMS-LIKE ALGORITHMS
%  Version    : v2.0
%  Author     : Zhengyi, Zhan
%  Release    : 9th Sep. 2021
%  Dependency : numCal.m (mLMS only)
%
%  USAGE
%  Creation of LMS
%  - lms = LMS(mode, fwhm, L, P);
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
%  - lms.SBR = SBR;
%  Set nVector
%  - lms.nVector = nVector;
%  Set P
%  - lms.P = P;
%
%  Method for LMS
%  - function [x, y] = lms.calLMS();
%    OUTPUT
%       x: estimated value of x
%       y: estimated value of y
%
%  Method for calculating beta for mLMS
%  - beta = lms.calBeta(k);
%    INPUT
%       k: the order of mLMS, in the range of [1,2]
%          the higher the value of k, the larger the unbiased est. area
%
%  Method for mLMS
%  - [x, y] = lms.calmLMS();
%    OUTPUT are as above
%
%  EXAMPLE
%
%  % Creation of an LMS object
%  lms = LMS(mode, fwhm, L, P);
%
%  % Set the nVector & SBR to be estimated before running algorithms
%  lms.nVector = n;
%  lms.SBR = SBR;
%
%  % Run LMS algorithm
%  [x, y] = lms.calLMS();
%
%  % Calculating beta for mLMS
%  beta = lms.calBeta(k);
%
%  % Run mLMS algorithm
%  [x, y] = lms.calmLMS();
%
%  Public Properties
%  - LMS.SBR
%  - LMS.nVector
%  - LMS.P

classdef LMS < handle
    properties (Access = private)
        mode
        fwhm
        L
        beta
    end
    
    properties (Access = public)
        P
        SBR
        nVector
    end
    
    methods
        
        function obj = LMS(mode, fwhm, L, P)
            
            narginchk(3,4);
            obj.mode = mode;
            obj.fwhm = fwhm;
            obj.L = L;
            
            if contains(mode, 'LDS')
                if length(P) > 3
                    error('The size of P must be smaller than 3.');
                else
                    obj.P = P;
                end
            end
            
        end
        
        function [x, y] = calLMS(obj)
            
            [x, y] = singleLMS(obj.mode, obj.nVector, obj.L, obj.fwhm, obj.P);
            
        end
        
        function betaF = calBeta(obj, k)
            
            N = sum(obj.nVector);
            [x,y] = meshgrid(-obj.L/2:obj.L/2,-obj.L/2:obj.L/2);
            M = numel(x);
            psf = numCal(obj.mode, obj.fwhm, x, y, obj.L);
            bg = sum(psf, 3) / (obj.SBR+1) / size(psf, 3);
            psfBg = psf + bg;
            lambda = psfBg ./ sum(psfBg, 3);
            p = lambda ./ sum(lambda, 3);
            
            if contains(obj.mode, 'doughnut', 'IgnoreCase', true)
                obj.P = 1;
            end
            
            [xLMS, yLMS, p0Sum] = multiLMS(obj.mode, p, obj.L, obj.fwhm, obj.P);
            
            options = optimoptions(@fminunc,'Display','off');
            switch k
                case 1
                    fun = @(beta) sum(sqrt(1/2/M*...
                        (x - xLMS .* (beta(1) +...
                        beta(2)*(N-1)/N*p0Sum)).^2 + ...
                        (y - yLMS .* (beta(1) +...
                        beta(2)*(N-1)/N*p0Sum)).^2),'all');
                    betaF = fminunc(fun,[0,0],options);
                case 2
                    fun = @(beta) sum(sqrt(1/2/M*...
                        (x - xLMS .* (beta(1) +...
                        (beta(2)+beta(3)/N)*(N-1)/N*p0Sum +...
                        beta(3)*(N-1)*(N-2)/N^2*p0Sum.^2)).^2 + ...
                        (y - yLMS .* (beta(1) +...
                        (beta(2)+beta(3)/N)*(N-1)/N*p0Sum +...
                        beta(3)*(N-1)*(N-2)/N^2*p0Sum.^2)).^2),'all');
                    betaF = fminunc(fun,[0,0,0],options);
                otherwise
                    error('The entered k value is not supported.');
            end
            
            obj.beta = betaF;
            
        end
        
        function [xmLMS, ymLMS] = calmLMS(obj)

            if isempty(obj.beta)
                error(['The calBeta function must be calculated before '...
                    'running calmLMS function.']);
            end
            
            [xLMS, yLMS, p0Sum] = singleLMS(obj.mode, obj.nVector, obj.L, obj.fwhm, obj.P);
            [xmLMS, ymLMS] = mLMS(obj.beta, sum(obj.nVector), p0Sum, xLMS, yLMS);
            
        end
        
    end
    
end

function [xLMS, yLMS, p0Sum] = singleLMS(mode, n, L, fwhm, P)

p = n / sum(n);
item = 1/(1-L^2*log(2)/fwhm^2);

if contains(mode, 'doughnut', 'IgnoreCase', true)
    if contains(mode, ["3", "4"])
        xTmp = -L/2*(-p(3)+(p(1)+p(2))/2);
        yTmp = -L/2*(sqrt(3)/2*(p(2)-p(1)));
        if contains(mode, '4')
            p0Sum = p(4);
        end
    elseif contains(obj.mode, ["6", "7"])
        xTmp = L/2*(-p(3)+p(6)+(p(1)-p(2)-p(4)+p(5))/2);
        yTmp = L/2*(sqrt(3)/2*(p(1)+p(2)-p(4)-p(5)));
        if contains(mode, '7')
            p0Sum = p(7);
        end
    else
        error('Mode error.');
    end
    
elseif contains(mode, 'LDS', 'IgnoreCase', true)
    Polar = PExtend(P, mode);
    nV = n ./ Polar;
    p = nV / sum(nV);
    
    if contains(mode, 'LDS-2', 'IgnoreCase', true)
        if contains(mode, ["5" "6"])
            xTmp = L/2*(p(5) - p(4));
            yTmp = L/2*(p(1) - p(2));
            if contains(mode, '6')
                p0Sum = p(3)+p(6);
            end
        elseif contains(mode, '4')
            xTmp = L/2*(p(4) - p(3));
            yTmp = L/2*(p(1) - p(2));
        else
            error('Mode error.');
        end
        
    elseif contains(mode, 'LDS-3', 'IgnoreCase', true)
        if contains(mode, '9')
            xTmp = sqrt(3)/4*L*(-p(4) + p(5) + p(7) - p(8));
            yTmp = L/2*(p(1) - p(2) + (p(4) - p(5) + p(7) - p(8))/2);
            p0Sum = p(3)+p(6)+p(9);
        elseif contains(mode, '7')
            xTmp = sqrt(3)/4*L*(-p(4) + p(5) + p(6) - p(7));
            yTmp = L/2*(p(1) - p(2) + (p(4) - p(5) + p(6) - p(7))/2);
        elseif contains(mode, '6')
            xTmp = sqrt(3)/4*L*(-p(3) + p(4) + p(5) - p(6));
            yTmp = L/2*(p(1) - p(2) + (p(3) - p(4) + p(5) - p(6))/2);
        else
            error('Mode error.');
        end
        
    else
        error('Mode error.');
    end
    
else
    error('Mode error.');
end

xLMS = item * xTmp;
yLMS = item * yTmp;

end

function [xLMS, yLMS, p0Sum] = multiLMS(mode, p, L, fwhm, P)

item = 1/(1-L^2*log(2)/fwhm^2);

if contains(mode, 'doughnut', 'IgnoreCase', true)
    if contains(mode, ["3", "6"])
        error(['Using calBeta function to calculate beta '...
            'for mLMS must have full centered exposes.']);
    elseif contains(mode, '4')
        xTmp = -L/2*(-p(:,:,3)+(p(:,:,1)+p(:,:,2))/2);
        yTmp = -L/2*(sqrt(3)/2*(p(:,:,2)-p(:,:,1)));
        p0Sum = p(:,:,4);
    elseif contains(mode, '7')
        xTmp = L/2*(-p(:,:,3)+p(:,:,6)+(p(:,:,1)-p(:,:,2)-p(:,:,4)+p(:,:,5))/2);
        yTmp = L/2*(sqrt(3)/2*(p(:,:,1)+p(:,:,2)-p(:,:,4)-p(:,:,5)));
        p0Sum = p(:,:,7);
    else
        error('Mode error.');
    end
    
elseif contains(mode, 'LDS', 'IgnoreCase', true)
    Polar = PExtend(P, mode);
    p = p ./ dotDimExpan(Polar);
    
    if contains(mode, 'LDS-2', 'IgnoreCase', true)
        if contains(mode, ["4" "5"])
            error(['Using calBeta function to calculate beta '...
                'for mLMS must have full centered exposes.']);
        elseif contains(mode, '6')
            xTmp = L/2*(p(:,:,5) - p(:,:,4));
            yTmp = L/2*(p(:,:,1) - p(:,:,2));
            p0Sum = p(:,:,3) + p(:,:,6);
        else
            error('Mode error.');
        end
        
    elseif contains(mode, 'LDS-3', 'IgnoreCase', true)
        if contains(mode, ["6" "7"])
            error(['Using calBeta function to calculate beta '...
                'for mLMS must have full centered exposes.']);
        elseif contains(mode, '9')
            xTmp = sqrt(3)/4*L*(-p(:,:,4) + p(:,:,5) + p(:,:,7) - p(:,:,8));
            yTmp = L/2*(p(:,:,1) - p(:,:,2) + (p(:,:,4) - p(:,:,5) + p(:,:,7) - p(:,:,8))/2);
            p0Sum = p(:,:,3) + p(:,:,6) + p(:,:,9);
        else
            error('Mode error.');
        end
        
    else
        error('Mode error.');
    end
    
else
    error('Mode error.');
end

xLMS = item * xTmp;
yLMS = item * yTmp;

end

function [xmLMS, ymLMS] = mLMS(beta, N, p0Sum, xLMS, yLMS)

switch length(beta)
    case 2
        xmLMS = xLMS .* (beta(1) + (N-1)/N*beta(2) * p0Sum);
        ymLMS = yLMS .* (beta(1) + (N-1)/N*beta(2) * p0Sum);
    case 3
        xmLMS = xLMS .* (beta(1) + (N-1)/N*(beta(2)+beta(3)/N) * p0Sum +...
            beta(3)*(N-1)*(N-2)/N^2 * p0Sum.^2);
        ymLMS = yLMS .* (beta(1) + (N-1)/N*(beta(2)+beta(3)/N) * p0Sum +...
            beta(3)*(N-1)*(N-2)/N^2 * p0Sum.^2);
end

end
