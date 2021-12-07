function out = normCustom(x, mode)
%NORMCUSTOM Custom normalization function
%   out = normCustom(x) normalize x to [0,1]
% -------------------------------------------------------------------------
% Coded by Zhengyi Zhan
% July 1, 2021
% -------------------------------------------------------------------------
narginchk(1,2);
if nargin < 2
    mode = 1;
end
switch mode
    case 1
        out = (x - min(x,[],'all')) / (max(x,[],'all') - min(x,[],'all'));
    case 2
        out = x / max(x,[],'all');
end
end

