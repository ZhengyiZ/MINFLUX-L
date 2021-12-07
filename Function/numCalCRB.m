function crb = numCalCRB(CRBFlag, mode, fwhm, L, N, SBR, x0, y0, A0, phi0)
%NUMCALCRB calculates CRB using numerical approximation
%   crb = numCalCRB(CRBFlag, mode, fwhm, L, N, SBR, x0, y0, A0, phi0)
% CRBFlag: 'xy' | 'x' | 'y' | 'A' | 'phi'
% A0 & alpha0 are only required in LDS modes
% -------------------------------------------------------------------------
% Coded by Zhengyi Zhan
% July 1, 2021
% -------------------------------------------------------------------------

narginchk(8,10);

[psf, psfGx, psfGy] = numCal(mode, fwhm, x0, y0, L);

switch CRBFlag
    case 'xy'
        if contains(mode, 'LDS', 'IgnoreCase', true)
            alpha = modeAlpha(mode);
            dotJ = cos(phi0 - alpha);
            if ndims(psf) == 3
                n = psf .* dotDimExpan(A0 .* dotJ.^2 + (1-A0)/2);
                dudx = dotDimExpan(A0 * dotJ.^2 + (1-A0)/2) .* psfGx;
                dudy = dotDimExpan(A0 * dotJ.^2 + (1-A0)/2) .* psfGy;
            else
                n = psf .* (A0 .* dotJ.^2 + (1-A0)/2);
                dudx = (A0 * dotJ.^2 + (1-A0)/2) .* psfGx;
                dudy = (A0 * dotJ.^2 + (1-A0)/2) .* psfGy;
            end
            u = n;
        else
            u = psf;
            dudx = psfGx;
            dudy = psfGy;
        end
        if ndims(u) == 3
            v = sum(u,3);
            if SBR == Inf
                p = u./v;
                dpdx = ( v.*dudx - u.*sum(dudx,3) ) ./ v.^2;
                dpdy = ( v.*dudy - u.*sum(dudy,3) ) ./ v.^2;
            else
                p = u./v * SBR/(SBR+1) + 1/size(psf,3)/(SBR+1);
                dpdx = ( v.*dudx - u.*sum(dudx,3) ) ./ v.^2 * SBR / (SBR+1);
                dpdy = ( v.*dudy - u.*sum(dudy,3) ) ./ v.^2 * SBR / (SBR+1);
            end
            crb = sqrt( 1 / (2*N) * sum(1./p .* (dpdx.^2+dpdy.^2),3) ./...
                ( sum(1./p .* dpdx.^2,3) .* sum(1./p .* dpdy.^2,3) - ...
                sum(1./p .* dpdx.*dpdy,3).^2) );
        else
            v = sum(u);
            if SBR == Inf
                p = u./v;
                dpdx = ( v.*dudx - u.*sum(dudx) ) ./ v.^2;
                dpdy = ( v.*dudy - u.*sum(dudy) ) ./ v.^2;
            else
                p = u./v * SBR/(SBR+1) + 1/length(psf)/(SBR+1);
                dpdx = ( v.*dudx - u.*sum(dudx) ) ./ v.^2 * SBR / (SBR+1);
                dpdy = ( v.*dudy - u.*sum(dudy) ) ./ v.^2 * SBR / (SBR+1);
            end
            crb = sqrt( 1 ./ (2*N) * sum(1./p .* (dpdx.^2+dpdy.^2)) ./...
                ( sum(1./p .* dpdx.^2) .* sum(1./p .* dpdy.^2) - ...
                sum(1./p .* dpdx.*dpdy).^2) );
        end
        
    case 'x'
        if contains(mode, 'LDS' , 'IgnoreCase', true)
            alpha = modeAlpha(mode);
            dotJ = cos(phi0 - alpha);
            n = psf .* (A0 .* dotJ.^2 + (1-A0)/2);
            u = n;
            du = (A0 * dotJ.^2 + (1-A0)/2) .* psfGx;
        else
            u = psf;
            du = psfGx;
        end
        v = sum(u);
        if SBR == Inf
            p = u./v;
            dp = ( v.*du - u.*sum(du) ) ./ v.^2;
        else
            p = u./v * SBR/(SBR+1) + 1/length(psf)/(SBR+1);
            dp = ( v.*du - u.*sum(du) ) ./ v.^2 * SBR / (SBR+1);
        end
        crb = sqrt( 1 ./ ( N * sum(1./p .* dp.^2) ) );
    case 'y'
        if contains(mode, 'LDS', 'IgnoreCase', true)
            alpha = modeAlpha(mode);
            dotJ = cos(phi0 - alpha);
            n = psf .* (A0 .* dotJ.^2 + (1-A0)/2);
            u = n;
            du = (A0 * dotJ.^2 + (1-A0)/2) .* psfGy;
        else
            u = psf;
            du = psfGy;
        end
        v = sum(u);
        if SBR == Inf
            p = u./v;
            dp = ( v.*du - u.*sum(du) ) ./ v.^2;
        else
            p = u./v * SBR/(SBR+1) + 1/length(psf)/(SBR+1);
            dp = ( v.*du - u.*sum(du) ) ./ v.^2 * SBR / (SBR+1);
        end
        crb = sqrt( 1 ./ ( N * sum(1./p .* dp.^2) ) );
    case 'A'
        if contains(mode,'LDS', 'IgnoreCase', true)
            alpha = modeAlpha(mode);
            dotJ = cos(phi0 - alpha);
            if ndims(psf) == 3
                n = psf .* dotDimExpan(A0 .* dotJ.^2 + (1-A0)/2);
                u = n;
                v = sum(u,3);
                du = psf/2 .* dotDimExpan(cos(2*(alpha-phi0)));
                if SBR == Inf
                p = u./v;
                dp = ( v.*du - u.*sum(du,3) ) ./ v.^2;
                else
                    p = u./v * SBR/(SBR+1) + 1/length(psf)/(SBR+1);
                    dp = ( v.*du - u.*sum(du,3) ) ./ v.^2 * SBR / (SBR+1);
                end
                crb = sqrt( 1 ./ ( N * sum(1./p .* dp.^2,3) ) );
            else
                n = psf .* (A0 .* dotJ.^2 + (1-A0)/2);
                u = n;
                v = sum(u);
                du = psf/2 .* cos(2*(alpha-phi0));
                if SBR == Inf
                    p = u./v;
                    dp = ( v.*du - u.*sum(du) ) ./ v.^2;
                else
                    p = u./v * SBR/(SBR+1) + 1/length(psf)/(SBR+1);
                    dp = ( v.*du - u.*sum(du) ) ./ v.^2 * SBR / (SBR+1);
                end
                crb = sqrt( 1 ./ ( N * sum(1./p .* dp.^2, 'omitnan') ) );
            end
        else
            error('Only LDS modes support A-CRB.');
        end
    case 'phi'
        if contains(mode,'LDS', 'IgnoreCase', true)
            alpha = modeAlpha(mode);
            dotJ = cos(phi0 - alpha);
            if ndims(psf) == 3
                n = psf .* dotDimExpan(A0 .* dotJ.^2 + (1-A0)/2);
                u = n;
                v = sum(u,3);
                du = A0 * psf .* dotDimExpan(sin(2*(alpha - phi0)));
                if SBR == Inf
                p = u./v;
                dp = ( v.*du - u.*sum(du,3) ) ./ v.^2;
                else
                    p = u./v * SBR/(SBR+1) + 1/length(psf)/(SBR+1);
                    dp = ( v.*du - u.*sum(du,3) ) ./ v.^2 * SBR / (SBR+1);
                end
                crb = sqrt( 1 ./ ( N * sum(1./p .* dp.^2,3) ) );
            else
                n = psf .* (A0 .* dotJ.^2 + (1-A0)/2);
                u = n;
                v = sum(u);
                du = A0 * psf .* sin(2*(alpha - phi0));
                if SBR == Inf
                    p = u./v;
                    dp = ( v.*du - u.*sum(du) ) ./ v.^2;
                else
                    p = u./v * SBR/(SBR+1) + 1/length(psf)/(SBR+1);
                    dp = ( v.*du - u.*sum(du) ) ./ v.^2 * SBR / (SBR+1);
                end
                crb = sqrt( 1 ./ ( N * sum(1./p .* dp.^2, 'omitnan') ) );
            end
        else
            error('Only LDS modes support phi-CRB.');
        end
    otherwise
        error('CRBFlag error.');
end

end
