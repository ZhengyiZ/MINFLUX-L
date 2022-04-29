function err = est2errV4angle(a, a0)
%EST2ERRV4ANGLE calculates error from estimation for angles (degree)
%   err = est2eb(a, a0)
% -------------------------------------------------------------------------
% Coded by Zhengyi Zhan
% Oct 25, 2021
% -------------------------------------------------------------------------

err = zeros(size(a));
% Correction for angles
for i = 1:length(a0)
    err(i) = a(i) - a0(i);
    if a0(i) >= 0 && abs(err(i)) >= 90
        err(i) = err(i) + 180;
    elseif a0(i) < 0 && abs(err(i)) >= 90
        err(i) = err(i) - 180;
    end
end