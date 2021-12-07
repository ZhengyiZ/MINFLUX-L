function cRes = fftconv3(x, m, shape, gpuFlag)
%FFTCONV3 Three dimensional convolution using FFT
%   cRes = FFTCONV3(x, m) performs the 3-D convolution of matrices X and M
%
%   cRes = FFTCONV3(x, m, shape) where SHAPE is a string returns a
%   subsection of the 3-D convolution with size specified by SHAPE:
%
%       'same'      - (default) returns the central part of the convolution
%                     that is the same size as X using zero padding
%
%       'full'      - returns the full 2-D convolution
% 
%       'valid'     - eturns only those parts of the convolution
%                     that can be computed without padding
%
%   cRes = FFTCONV3(x, m, shape, gpuFlag) where gpuFlag is an integer 
%   controling whether to use GPU
%        1 for using GPU
%        0 for not using GPU
%
% -------------------------------------------------------------------------
% Coded by Liu, Xin
% Jun 10, 2021
% Updated by Zhan, Zhengyi
% Jul 12, 2021
% -------------------------------------------------------------------------

narginchk(2,4);
if nargin < 3
    shape = 'same';
    gpuFlag = 0;
elseif nargin < 4
    if isnumeric(shape)
        gpuFlag = shape;
        shape = 'same';
    else
        gpuFlag = 0;
    end
end
if ~ismember(shape, {'full' 'same' 'valid'})
    error("Shape must be one of 'full', 'same' and 'valid'.");
end
    
% matrix expansion
[x1, x2, x3] = size(x);
[m1, m2, m3] = size(m);

c = zeros(x1+m1-1, x2+m2-1, x3+m3-1);

% gpu available memory detection
if gpuFlag == 1
    am = gpuDevice(1).AvailableMemory;
    if am < 0.5*4*numel(c)
        fprintf('Available GPU memory < 50 percent of required, GPU computing of FFT has been turned off.');
        gpuFlag = 0;
    end
end

if gpuFlag == 1
    c = gpuArray(c);
end

d = c;

c(1:x1, 1:x2, 1:x3) = x;
C = fftn(c);
clear c;

d(1:m1, 1:m2, 1:m3) = m;
D = fftn(d);
clear d;

E = C.*D;
clear C D;

E = abs(ifftn(E));

switch shape
    case 'same'
        cRes = E(floor(m1/2)+1:floor(m1/2)+x1,...
            floor(m2/2)+1:floor(m2/2)+x2,...
            floor(m3/2)+1:floor(m3/2)+x3);
    case 'full'
        cRes = E;
    case 'valid'
        if x1 < m1 || x2 < m2 || x3 < m3
            cRes = [];
        else
            cRes = E(m1:x1, m2:x2, m3:x3);
        end
end

if gpuFlag
    cRes = gather(cRes);
end
