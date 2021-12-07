function alpha = modeAlpha(mode)
%MODEALPHA returns the alpha vector by mode
%   alpha = modeAlpha(mode)
%
%   INPUT PARAMETERS
%      mode: 'LDS-2-4' | 'LDS-2-5' | 'LDS-2-6'
%            'LDS-3-6' | 'LDS-3-7' | 'LDS-3-9'
% -------------------------------------------------------------------------
% Coded by Zhengyi Zhan
% Jun 29, 2021
% -------------------------------------------------------------------------
if contains(mode, ["gaussian", "doughnut"], 'IgnoreCase', true)
    error('modeTheta only support LDS modes.');
end

if contains(mode, 'LDS-2', 'IgnoreCase', true)
    
    if contains(mode, '4')
        alpha(1:2) = deg2rad(0);
        alpha(3:4) = deg2rad(90);
    elseif contains(mode, '5')
        alpha(1:3) = deg2rad(0);
        alpha(4:5) = deg2rad(90);
    elseif contains(mode, '6')
        alpha(1:3) = deg2rad(0);
        alpha(4:6) = deg2rad(90);
    else
        error('mode error.');
    end
    
elseif contains(mode, 'LDS-3', 'IgnoreCase', true)
    
    if contains(mode, '6')
        alpha(1:2) = deg2rad(0);
        alpha(3:4) = deg2rad(60);
        alpha(5:6) = deg2rad(-60);
    elseif contains(mode, '7')
        alpha(1:3) = deg2rad(0);
        alpha(4:5) = deg2rad(60);
        alpha(6:7) = deg2rad(-60);
    elseif contains(mode, '9')
        alpha(1:3) = deg2rad(0);
        alpha(4:6) = deg2rad(60);
        alpha(7:9) = deg2rad(-60);
    else
        error('mode error.');
    end
    
else
    
    error('mode error.');
    
end
alpha = alpha';
end
