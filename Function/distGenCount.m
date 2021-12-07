function n = distGenCount(pDist, nTotal)
%%DISTGENCOUNT generates photon counts based on probability distribution
% n = distGenCount(pDist, nTotal)
% -------------------------------------------------------------------------
% Coded by Zhengyi Zhan
% Sept 9, 2021
% -------------------------------------------------------------------------

% generate random number based on probability distribution
rs = random(pDist, nTotal, 1);

% count every channel
tb = tabulate(rs);
n = tb(:,2);

% in case the last few channels do not have even one count
while length(n) ~= length(pDist.Probabilities)
    n(end+1) = 0;
end

if isrow(n)
    n = n';
end

end
