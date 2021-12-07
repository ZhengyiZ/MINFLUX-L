function [err, bias] = est2eb(xV, x0, yV, y0)
%EST2EB calculates error & bias from estimation
%   [err, bias] = est2eb(xV, x0, yV, y0)
%   If xV (yV) is not a vector, this function will calculate along the last
%   dimension
% -------------------------------------------------------------------------
% Coded by Zhengyi Zhan
% Sept 9, 2021
% -------------------------------------------------------------------------

if nargin == 2
    if isvector(xV)
        err = sqrt( mean((xV-x0).^2) );
        bias = sqrt( (mean(xV)-x0).^2 );
    else
        err = sqrt( mean((xV-x0).^2, ndims(xV)) );
        bias = sqrt( (mean(xV, ndims(xV))-x0).^2 );
    end
elseif nargin == 4
    if isvector(xV)
        err = sqrt( 1/2 * mean((xV-x0).^2 + (yV-y0).^2) );
        bias = sqrt( 1/2 * ((mean(xV)-x0).^2 + (mean(yV)-y0).^2));
    else
        err = sqrt( 1/2 * mean((xV-x0).^2 + (yV-y0).^2, ndims(xV)) );
        bias = sqrt( 1/2 * ((mean(xV, ndims(xV))-x0).^2 + (mean(yV, ndims(yV))-y0).^2));
    end
else
    error('Not support more than two variables.');
end

end