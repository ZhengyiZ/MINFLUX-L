function out = dotDimExpan4(in, iLen, jLen, kLen)
%DOTDIMEXPAN4 Expan input array dimension
% out = dotDimExpan4(in) expans the array to iLen x jLen x kLen x 3
% -------------------------------------------------------------------------
% Coded by Zhengyi Zhan
% July 1, 2021
% -------------------------------------------------------------------------

out = zeros(iLen, jLen, kLen, 3);
for i = 1:iLen
    for j = 1:jLen
        for k = 1:kLen
            out(i,j,k,:) = in;
        end
    end
end

end
