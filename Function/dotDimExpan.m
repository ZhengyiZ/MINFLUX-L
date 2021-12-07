function out = dotDimExpan(in, len)
%DOTDIMEXPAN Expan input array dimension
% out = dotDimExpan(in) expans the array to 1 x 1 x n
% -------------------------------------------------------------------------
% Coded by Zhengyi Zhan
% July 1, 2021
% -------------------------------------------------------------------------
narginchk(1,2);

if length(in) ~= numel(in)
    error('The input argument must be an 1-D array.')
end

if nargin < 2
    out = zeros(1,1,length(in));
    for i = 1:length(in)
        out(:,:,i) = in(i);
    end
else
    out = zeros(len,3);
    for i = 1:len
        out(i,:) = in;
    end
end

end
