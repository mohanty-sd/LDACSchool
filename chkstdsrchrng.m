function ptFlags = chkstdsrchrng(xVec)
%Checks for points that are outside the standardized search range 
%V = CHKSTDSRCHRNG(X)
%Returns an array of logical indices V corresponding to valid/invalid
%points in X. A row (point) of X is invalid if any of the coordinates
%(columns for that row) fall outside the closed interval [0,1].
%Do Y = X(V,:) to retrieve only the valid rows or Y = X(~V,:) to retrieve
%invalid rows.

%Soumya D. Mohanty
%April 2012

[nrows,ncols]=size(xVec);
validPts = ones(1,nrows);
for lp = 1:nrows
    x = xVec(lp,:);
    if any(x<0|x>1)
        %Mark point as invalid
        validPts(lp) = 0;
    end
end
ptFlags = logical(validPts);
