function dispMultiInfo(mode, SBR, nTotal, poissFlag, methodFlag,...
    lr, lrPenalty, iterations, Amin, Amax, thetaMin, thetaMax, alphaFixed)

narginchk(6,14);

fprintf(['Under ' mode ' mode. ']);
fprintf('SBR = %d.\n', SBR);

if contains(mode, ["1D-2" "1D-3"])
    if Amin > Amax
        error('The Amin must be smaller than Amax.');
    elseif Amin == Amax
        fprintf('Fixed A = %d. ', Amin);
    else
        fprintf('Random A, range from %.2f to %.2f. \n', Amin, Amax);
    end
    if thetaMin > thetaMax
        error('The thetaMin must be smaller than thetaMax.');
    elseif thetaMin == thetaMax
        fprintf('Fixed theta = %d°. ', thetaMin);
    else
        fprintf('Random theta, range from %d° to %d°. \n', thetaMin, thetaMax);
    end
    if exist('alphaFixed', 'var')
        fprintf('Fixed α = %d. \n', alphaFixed);
    end
end

if poissFlag
    fprintf('Photon counts obey Poisson distribution ');
else
    fprintf('Photon counts directly take the means of Poisson distribution ');
end

fprintf('and N = %d. \n', nTotal);

switch methodFlag
    case 1
        fprintf('Using lr = %.3f, lrPenalty = %d, iterations = %d as Adam parameters.\n',...
            lr, lrPenalty, iterations);
    case 2
        
end

end
