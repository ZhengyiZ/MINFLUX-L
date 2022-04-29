function dispPSF(psf, mode, x0, y0)
%DISPPSF displays PSF by mode
% dispPSF(psf, mode, x0, y0)
% x0 & y0 need to be vectors, and are not necessary
% -------------------------------------------------------------------------
% Coded by Zhengyi Zhan
% Jun 9, 2021
% -------------------------------------------------------------------------
narginchk(1,4);
if nargin < 4
    x0 = 1:size(psf,2);
    y0 = 1:size(psf,1);
    if nargin < 2
        mode = 'default';
    end
end

if ndims(psf) ~= 3
    error('PSF must be a 3D martix.');
end

if ~isvector(x0) || ~isvector(y0)
    error('x0 & y0 must be vectors.');
end

psf = normCustom(psf);

top_margin = 0.06; % top margin
btm_margin = 0.06; % bottom margin
left_margin = 0.06;% left margin
right_margin = 0.1;% right margin
fig_margin = 0.08; % margin beween figures(sub)
clim = [min(psf,[],'all') max(psf,[],'all')];

if contains(mode, 'LDS')
    flag = 0;
    switch size(psf, 3)
        case 4 % LDS-2-4
            row = 2;
            col = 2;
            skip = [];
            seq = [1 2 3 4];
        case 5 % LDS-2-5
            row = 2;
            col = 3;
            skip = 6;
            seq = [1 2 3 4 5];
        case 6 % LDS-2-6 | LDS-3-6 | Doughnut-6
            row = 2;
            col = 3;
            skip = [];
            if contains(mode, '-3')
                seq = [1 3 5 2 4 6];
            else
                seq = [1 2 3 4 5 6];
            end
        case 7 % LDS-3-7
            row = 3;
            col = 3;
            skip = [6 9];
            seq = [1 2 3 4 5 6 7];
        case 9 % LDS-3-9
            row = 3;
            col = 3;
            skip = [];
            seq = 1:9;
        otherwise
            error('The size of PSF is not supported.');
    end
else
    row = 3;
    col = 3;
    flag = 1;
    switch size(psf,3)
        case 3
            skip = [1 3 5 6 7 9];
            seq = [1 3 2];
        case 4
            skip = [1 3 6 7 9];
            seq = [1 3 4 2];
        case 6
            skip = [3 5 9];
            seq = [1 2 6 3 5 4];
        case {7,14}
            skip = [3 9];
            seq = [1 2 6 7 3 5 4];
    end
    
end

count = 1;
% Calculate figure height and width according to rows and cols
fig_h = (1- top_margin - btm_margin - (row-1) * fig_margin) / row;
fig_w = (1 - left_margin - right_margin - (col-1) * fig_margin) / col;
if flag
    fig_s = (1 - left_margin - right_margin - fig_w*2 - fig_margin) / 2;
else
    fig_s = 0;
end
figure,
for i = 1 : row
    for j = 1 : col
        if ~sum(skip == (i-1)*col+j)
            % figure position: you can refer to 'help axes' to review the
            % parameter meaning, note that original point is lower-left
            if i == 1 || i == 3
                position = [left_margin + fig_s + (j-1)*(fig_margin+fig_w), ...
                    1- (top_margin + i * fig_h + (i-1) * fig_margin), ...
                    fig_w, fig_h];
            else
                position = [left_margin + (j-1)*(fig_margin+fig_w), ...
                    1- (top_margin + i * fig_h + (i-1) * fig_margin), ...
                    fig_w, fig_h];
            end
            axes('position', position);
            % draw colorful pictures...
            imagesc(x0,y0,psf(:,:,seq(count)), clim), axis image; axis off;
            count = count + 1;
        end
    end
end

% suptitle('EBP');

% draw colorbar
if contains(mode, ["doughnut-3" "gaussian-3" "doughnut-4" "gaussian-4"], 'IgnoreCase', true)
    cbPos = [1-right_margin-fig_margin-fig_s, btm_margin, 0.2, 1-(top_margin+btm_margin)];
else
    cbPos = [1-right_margin-fig_margin, btm_margin, 0.2, 1-(top_margin+btm_margin)];
end
axes('position', cbPos);
axis off;
colorbar();caxis(clim);
colormap('hot');

end
