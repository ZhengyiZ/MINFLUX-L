function [err, bias] = est2eb4angle(a, a0)
%EST2EB4ANGLE calculates error & bias from estimation for angles (degree)
%   [err, bias] = est2eb(a, a0)
%   If xV (yV) is not a vector, this function will calculate along the last
%   dimension
% -------------------------------------------------------------------------
% Coded by Zhengyi Zhan
% Sept 9, 2021
% -------------------------------------------------------------------------

% Correction for angles
if a0 >= 0
    a(abs(a-a0)>=90) = a(abs(a-a0)>=90) + 180;
elseif a0 < 0
    a(abs(a-a0)>=90) = a(abs(a-a0)>=90) - 180;
end

if isvector(a)
    err = sqrt( mean((a-a0).^2) );
    bias = sqrt( (mean(a)-a0).^2 );
else
    err = sqrt( mean((a-a0).^2, ndims(a)) );
    bias = sqrt( (mean(a, ndims(a))-a0).^2 );
end

end