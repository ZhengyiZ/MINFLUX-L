function [n, SBRest, count] = phyGenCount(p0, bg, N, laserPower)
%%PHYGENCOUNT generates photon counts based on physical process
% [n, SBRest, count] = phyGenCount(p0, bg, N, laserPower)
% INPUT PARAMETERS
% p0: the emission probabilities ignoring dark count
% bg: the backgroud
% N: the limit of total number of photons
% laserPower: a constant controls the number of summation
% OUTPUT PARAMETERS
% n: photons collected each exposure, vector
% SBRest: the estimated SBR from the residual background random numbers
% count: the number of EBP used for collection, scalar
% -------------------------------------------------------------------------
% Coded by Zhengyi Zhan
% July 10, 2021
% -------------------------------------------------------------------------
if ndims(p0) == 3
    p0 = squeeze(p0);
end
if isrow(p0)
    p0 = p0';
end

p0 = p0 * laserPower/10*(log(N)+1);
bg = bg * laserPower/10*(log(N)+1);

pData = zeros(length(p0), 1000*length(p0)/laserPower);
bgData = pData;
for i = 1:length(p0)
    pData(i,:) = poissrnd(p0(i), 1000*length(p0)/laserPower, 1);
    bgData(i,:) = poissrnd(bg, 1000*length(p0)/laserPower, 1);
end
pb = pData + bgData;
n = zeros(size(p0));
count = 1;
while sum(n) < N
    if sum(n+pb(:,count)) > 1.2*N && sum(n) > 0
        count = count - 1;
        break;
    else
        n = n + pb(:, count);
    end
    count = count + 1;
    if count > length(pData)
        error('not enough data.');
    end
end

SBRest = sum(n)/mean(sum(bgData(:,count+1:end),1))/count - 1;

end
