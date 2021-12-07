function dispInfo(mode, SBR, nTotal, methodFlag,...
    lr, lrPenalty, iterations)

narginchk(4,7);

fprintf(['Under ' mode ' mode. ']);
fprintf('SBR = %d ', SBR);
fprintf('and N = %d. \n', nTotal);

switch methodFlag
    case 1
        fprintf('Using lr = %.3f, lrPenalty = %d & iterations = %d as Adam parameters.\n',...
            lr, lrPenalty, iterations);
    case 2
        
end

end
