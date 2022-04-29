function n = phyGenCount1(p0, bg, laserPower)
%%PHYGENCOUNT generates poissonian photon counts only once
% n = phyGenCount(p0, bg, laserPower)
% n: photons collected each channel, vector
% laserPower: <nTotal> ~ laserPower/10
% -------------------------------------------------------------------------
% Coded by Zhengyi Zhan
% Oct 25, 2021
% -------------------------------------------------------------------------
if ndims(p0) == 3
    p0 = squeeze(p0);
end
if isrow(p0)
    p0 = p0';
end

p = (p0 + bg) / sum(p0+bg);
p = p * laserPower / 10;
n = p;

for i = 1:length(p)
    n(i) = poissrnd(p(i), 1, 1);
end

end
