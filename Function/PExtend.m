function Pout = PExtend(Pin, mode)
%PEXTEND Extend the length of P vector to the length required by mode
% Pout = PExtend(Pin, mode) extend the length of Pin by mode
% -------------------------------------------------------------------------
% Coded by Zhengyi Zhan
% July 1, 2021
% -------------------------------------------------------------------------

if contains(mode, 'LDS-2')
    if length(Pin) < 2
        error('The size of P vector must be 2 in LDS-2 modes.');
    elseif length(Pin) > 2
        Pin = Pin(1:2);
    end
    switch mode
        case 'LDS-2-4'
            Pout = zeros(4,1);
            Pout(1:2) = Pin(1);
            Pout(3:4) = Pin(2);
        case 'LDS-2-5'
            Pout = zeros(5,1);
            Pout(1:3) = Pin(1);
            Pout(4:5) = Pin(2);
        case 'LDS-2-6'
            Pout = zeros(6,1);
            Pout(1:3) = Pin(1);
            Pout(4:6) = Pin(2);
    end
elseif contains(mode, 'LDS-3')
    if length(Pin) < 3
        error('The size of P vector must be 3 in LDS-3 modes.');
    elseif length(Pin) > 3
        Pin = Pin(1:3);
    end
    switch mode
        case 'LDS-3-6'
            Pout = zeros(6,1);
            Pout(1:2) = Pin(1);
            Pout(3:4) = Pin(2);
            Pout(5:6) = Pin(3);
        case 'LDS-3-7'
            Pout = zeros(7,1);
            Pout(1:3) = Pin(1);
            Pout(4:5) = Pin(2);
            Pout(6:7) = Pin(3);
        case 'LDS-3-9'
            Pout = zeros(9,1);
            Pout(1:3) = Pin(1);
            Pout(4:6) = Pin(2);
            Pout(7:9) = Pin(3);
    end
end

