function [x, y, orient] = genLSC(a, b, p, q, phi, t)
%GENLSC Lissajous curve at every t moment
%   lsc returns the x and y value of Lissajous curve
%   x = a*sin(p*t+phi), y = b*sin(q*t)
%   [x,y,orient] = lsc(a, b, p, q, phi, t)
% -------------------------------------------------------------------------
% Coded by Zhengyi Zhan
% July 1, 2021
% -------------------------------------------------------------------------
x = a * sin( p*t + phi );
y = b * sin( q*t );
dx = a * p * cos( p*t + phi );
dy = b * q * cos( q*t );
orient = atan(dy./dx);

end
