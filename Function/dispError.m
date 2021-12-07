function dispError(methodFlag, evalFlag, mode, x0, xEst, y0, yEst,...
    A0, AEst, theta, alpha, alphaEst, n, fwhm, L, SBR, P)

narginchk(7,17);

if evalFlag
    narginchk(11,17);
    if nargin == 11
        n = A0;
        fwhm = AEst;
        L = theta;
        SBR = alpha;
    end
end

switch methodFlag
    case 1
        fprintf('Adam: \n');
    case 2
        fprintf('Grid Search： \n');
    case 3
        fprintf('LMS： \n');
    case 4
        fprintf('mLMS： \n');
    case 5
        fprintf('numLMS： \n');
    otherwise
        erorr('methodFlag error.');
end

switch methodFlag
    case {1, 2}
        if contains(mode, ["gaussian", "doughnut"])
            fprintf('x0 = %6.2f,\t y0 = %6.2f\n', x0, y0);
            fprintf('xE = %6.2f,\t yE = %6.2f,\t RMSE  = %6.2f\n',...
                abs(x0-xEst), abs(y0-yEst),...
                sqrt(1/2*((x0-xEst)^2+(y0-yEst)^2)));
            if evalFlag
                fprintf('evalFun = %6.2f\n',...
                    evalFun(mode,n,numCal(mode,fwhm,xEst,yEst,L), SBR));
            end
        else
            alphaEst = rad2deg(alphaEst);
            fprintf('x0 = %6.2f,\t y0 = %6.2f,\t A  = %5.2f,\t α  = %6.2f,\t theta  = %6.2f,\n',...
                x0, y0, A0, alpha, theta);
            if alphaEst - alpha > 90
                alphaEst = alphaEst - 180;
            elseif alphaEst - alpha < -90
                alphaEst = alphaEst + 180;
            end
            fprintf('xE = %6.2f,\t yE = %6.2f,\t AE = %5.2f,\t αE = %6.2f,\t RMSE   = %6.2f\n',...
                abs(x0-xEst), abs(y0-yEst), abs(A0-AEst),...
                abs(alpha-alphaEst), sqrt(1/2*((x0-xEst)^2+(y0-yEst)^2)));
            if evalFlag
                fprintf('evalFun = %6.2f\n',...
                    evalFun(mode, n, numCal(mode,fwhm,xEst,yEst,L), SBR, ...
                    PExtend(P,mode), AEst, alphaEst));
            end
        end
    case {3, 4, 5}
        if contains(mode, ["gaussian", "doughnut"])
            fprintf('x0 = %6.2f,\t y0 = %6.2f\n', x0, y0);
            fprintf('xE = %6.2f,\t yE = %6.2f,\t RMSE  = %6.2f\n',...
                abs(x0-xEst), abs(y0-yEst),...
                sqrt(1/2*((x0-xEst)^2+(y0-yEst)^2)));
            if evalFlag
                fprintf('evalFun = %6.2f\n',...
                    evalFun(mode,n,numCal(mode,fwhm,xEst,yEst,L), SBR));
            end
        else
            fprintf('x0 = %6.2f,\t y0 = %6.2f\n', x0, y0);
            fprintf('xE = %6.2f,\t yE = %6.2f,\t RMSE  = %6.2f\n',...
                abs(x0-xEst), abs(y0-yEst),...
                sqrt(1/2*((x0-xEst)^2+(y0-yEst)^2)));
        end
    otherwise
        error('methodFlag error.');
end

end
