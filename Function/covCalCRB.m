function [crbXY, crbA, crbPhi] = covCalCRB(CovMat)

crbXY = sqrt((CovMat(1,1)+CovMat(2,2))/2);
crbA = sqrt(CovMat(3,3));
crbPhi = rad2deg(sqrt(CovMat(4,4)));

end
