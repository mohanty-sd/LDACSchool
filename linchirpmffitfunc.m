function [fitVal,varargout] = linchirpmffitfunc(xVec,params)
%Fitness function for matched filtering a linear chirp
%F = LINCHIRPMFFITFUNC(X,P)
%Compute the matched filter output maximum for the linear chirp
%parameters in each row.  The fitness values are returned in F.
%X is standardized, that is 0<=X(i,j)<=1, and the array of minimum and
%maximum ('rmin' and 'rmax' respectively) values in P are used to convert X(i,j)
%internally before computing fitness: 
%X(:,j) -> X(:,j)*(rmax(j)-rmin(j))+rmin(j).
%The other fields of P are as follows:
% 'fftData' : FFT of the data vector
% 'timeVec' : The sampling times
% 'fs' : The sampling frequency
% 'sigLen' : The length of the chirp signal (and templates)
%
%For standardized coordinates, F = infty if a point X(i,:) falls
%outside the hypercube defined by 0<=X(i,j)<=1.
%
%[F,R] =  LINCHIRPMFFITFUNC(X,P)
%returns the real coordinates in R. 
%
%[F,R,Xp] = LINCHIRPMFFITFUNC(X,P)
%Returns the standardized coordinates in Xp. This option is to be used when
%there are special boundary conditions (such as wrapping of angular
%coordinates) that are better handled by the fitness function itself.

%Soumya D. Mohanty, Aug 2015
%Just a renamed version of the rastrigin benchmark function.

%Soumya D. Mohanty
%June, 2011
%April 2012: Modified to switch between standardized and real coordinates.

%Shihan Weerathunga
%April 2012: Modified to add the function rastrigin.

%Soumya D. Mohanty
%May 2016: New optional output argument introduced in connection with
%handling of special boundary conditions.

%Soumya D. Mohanty
%Dec 2017: Modified LDACPSOTESTFUNC.
%==========================================================================

%rows: points
%columns: coordinates of a point
[nrows,~]=size(xVec);

%storage for fitness values
fitVal = zeros(nrows,1);
validPts = ones(nrows,1);

%Check for out of bound coordinates and flag them
validPts = chkstdsrchrng(xVec);
%Set fitness for invalid points to infty
fitVal(~validPts)=inf;
xVec(validPts,:) = s2rvector(xVec(validPts,:),params);

%-------------Quantities specified to this fitness function--------
nSamples = length(params.timeVec);
templtNSamples = floor(params.sigLen*params.fs);
q0Vec = zeros(1,nSamples);
q1Vec = zeros(1,nSamples);
templateTimeVec = params.timeVec(1:templtNSamples);
%------------------------------------------------------------------

for lpc = 1:nrows
    if validPts(lpc)
    % Only the body of this block should be replaced for different fitness
    % functions
        x = xVec(lpc,:);
        %Compute the negative of the matched filter output (because the
        %LDACPSO function is a minimizer) and return its minimum
        q0Vec(:)=0;
        q1Vec(:)=0;
        q0Vec(1:templtNSamples) = sin(2*pi*(x(1)*templateTimeVec+...
                                   x(2)*templateTimeVec.^2));
        q0Vec = q0Vec/norm(q0Vec);
        q1Vec(1:templtNSamples) = cos(2*pi*(x(1)*templateTimeVec+...
            x(2)*templateTimeVec.^2));
        q1Vec = q1Vec/norm(q1Vec);
        % Fourier transforms of the templates
        fftq0 = fft(q0Vec);
        fftq1 = fft(q1Vec);
        % We do the correlation in the Fourier domain.
        q0Corr = real(ifft(params.fftData.*conj(fftq0)));
        q1Corr = real(ifft(params.fftData.*conj(fftq1)));
        %% Sum of squares of quadrature correlation outputs
        mfOut = sqrt(q0Corr.^2+q1Corr.^2);
        fitVal(lpc) = -max(mfOut);
    end
end

%Return real coordinates if requested
if nargout > 1
    varargout{1}=xVec;
    if nargout > 2
            varargout{2} = r2svector(xVec,params);
    end
end

function rVec = s2rvector(xVec,params)
%Convert standardized coordinates to real using non-uniform range limits
%R = S2RSCALAR(X,P)
%Takes standardized coordinates in X (coordinates of one point per row) and
%returns real coordinates in R using the range limits defined in P.rmin and
%P.rmax. The range limits can be different for different dimensions. (If
%they are same for all dimensions, use S2RSCALAR instead.)

%Soumya D. Mohanty
%April 2012

[nrows,ncols] = size(xVec);
rVec = zeros(nrows,ncols);
rmin = params.rmin;
rmax = params.rmax;
rngVec = rmax-rmin;
for lp = 1:nrows
    rVec(lp,:) = xVec(lp,:).*rngVec+rmin;
end

function xVec = r2svector(rVec,params)
%Convert real coordinates to standardized ones.
%X = R2SVECTOR(R,P)
%Takes standardized coordinates in X (coordinates of one point per row) and
%returns real coordinates in R using the range limits defined in P.rmin and
%P.rmax. The range limits can be different for different dimensions. (If
%they are same for all dimensions, use S2RSCALAR instead.)

%Soumya D. Mohanty
%May 2016
[nrows,ncols] = size(rVec);
xVec = zeros(nrows,ncols);
rmin = params.rmin;
rmax = params.rmax;
rngVec = rmax-rmin;
for lp = 1:nrows
    xVec(lp,:) = (rVec(lp,:)-rmin)./rngVec;
end

function rVec = s2rscalar(xVec,rngMin,rngMax)
%Convert standardized coordinates to real using uniform range limits
%R = S2RSCALAR(X,R1,R2)
%Takes standardized coordinates in X (coordinates of one point per row) and
%returns real coordinates in R using the range limits R1 < R2.

%Soumya D. Mohanty
%April 2012

%Number of rows = number of points
%Number of columns = number of dimensions
[nrows,ncols]=size(xVec);
%Storage for real coordinates
rVec = zeros(nrows,ncols);
%Range for each coordinate dimension is the same
rmin = rngMin*ones(1,ncols);
rmax = rngMax*ones(1,ncols);
rngVec = rmax-rmin;
%Apply range to standardized coordinates and convert to real coordinates
for lp = 1:nrows
    rVec(lp,:) = xVec(lp,:).*rngVec+rmin;
end

function xVec = r2sscalar(rVec,rngMin,rngMax)
%Convert real coordinates to standardized ones using uniform range limits
%R = R2SSCALAR(X,R1,R2)
%Takes standardized coordinates in X (coordinates of one point per row) and
%returns real coordinates in R using the range limits R1 < R2.

%Soumya D. Mohanty
%May 2016

%Number of rows = number of points
%Number of columns = number of dimensions
[nrows,ncols]=size(rVec);
%Storage for real coordinates
xVec = zeros(nrows,ncols);
%Range for each coordinate dimension is the same
rmin = rngMin*ones(1,ncols);
rmax = rngMax*ones(1,ncols);
rngVec = rmax-rmin;
%Apply range to standardized coordinates and convert to real coordinates
for lp = 1:nrows
    xVec(lp,:) = (rVec(lp,:)-rmin)./rngVec;
end