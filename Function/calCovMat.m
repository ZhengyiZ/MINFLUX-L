function CovMat = calCovMat(psf,psfGx,psfGy,A0,alpha,phi0,SBR,N)

absEff = (A0 * cosd(alpha - phi0).^2 + (1-A0)/2);

if ndims(psf) > 1
    psf = squeeze(psf);
    psfGx = squeeze(psfGx);
    psfGy = squeeze(psfGy);
end

u = psf .* absEff;
v = sum(u);

dudx = absEff .* psfGx;
dudy = absEff .* psfGy;
dudA = psf/2 .* cosd(2*(alpha-phi0));
dudPhi = A0 * psf .* sind(2*(alpha-phi0));

if SBR == Inf
    p = u./v;
    dpdx = ( v.*dudx - u.*sum(dudx) ) ./ v.^2;
    dpdy = ( v.*dudy - u.*sum(dudy) ) ./ v.^2;
    dpdA = ( v.*dudA - u.*sum(dudA) ) ./ v.^2;
    dpdPhi = ( v.*dudPhi - u.*sum(dudPhi) ) ./ v.^2;
else
    p = SBR/(SBR+1)*u./v + 1/(SBR+1)/9;
    dpdx = ( v.*dudx - u.*sum(dudx) ) ./ v.^2 * SBR / (SBR+1);
    dpdy = ( v.*dudy - u.*sum(dudy) ) ./ v.^2 * SBR / (SBR+1);
    dpdA = ( v.*dudA - u.*sum(dudA) ) ./ v.^2 * SBR / (SBR+1);
    dpdPhi = ( v.*dudPhi - u.*sum(dudPhi) ) ./ v.^2 * SBR / (SBR+1);
end
F = [sum(dpdx.^2./p)        sum(dpdx.*dpdy./p)      sum(dpdx.*dpdA./p)      sum(dpdx.*dpdPhi./p);
     sum(dpdx.*dpdy./p)     sum(dpdy.^2./p)         sum(dpdy.*dpdA./p)      sum(dpdy.*dpdPhi./p);
     sum(dpdx.*dpdA./p)     sum(dpdy.*dpdA./p)      sum(dpdA.^2./p)         sum(dpdA.*dpdPhi./p);
     sum(dpdx.*dpdPhi./p)   sum(dpdy.*dpdPhi./p)    sum(dpdA.*dpdPhi./p)    sum(dpdPhi.^2./p)];
F = F * N;
CovMat = inv(F);

end
