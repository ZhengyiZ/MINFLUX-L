function cmap = genHSVCmap(h, vl)
%GENHSVCMAP generates cmap based on HSV
%   cmap = genHSVCmap(h, vl)
%   vl: the highest value of lightness
% -------------------------------------------------------------------------
% Coded by Zhengyi Zhan
% July 1, 2021
% -------------------------------------------------------------------------
hsvcmap = [ones(256,1)*h ones(256,1) ones(256,1)];
for i = 1:256
    hsvcmap(i,2) = (i-1)/255;
    hsvcmap(i,3) = vl+(256-i)/255*(1-vl);
end
cmap = hsv2rgb(hsvcmap);
end