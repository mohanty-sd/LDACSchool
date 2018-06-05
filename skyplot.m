function [] = skyplot(alphaVec,deltaVec,fHandle)
%Plot a function of sky angles on the unit sphere
%SKYPLOT(A,D,F)
%A and D are the vector of azimuthal and polar angle values respectively. F
%is the handle to the function of A and D that is to be plotted on the
%sphere. The polar angle range is 0 (z-axis) to 180 (negative z-axis). The
%angles are assumed to be in radian.

%Soumya D. Mohanty, June 2018

%Generate the X, Y, Z meshes corresponding to the azimuthal and polar
%angles
[A,D] = meshgrid(alphaVec,deltaVec);
X = sin(D).*cos(A);
Y = sin(D).*sin(A);
Z = cos(D);

%Generate function values
fVals = zeros(length(deltaVec),length(alphaVec));
for lp1 = 1:length(alphaVec)
    for lp2 = 1:length(deltaVec)
        fVals(lp2,lp1) = fHandle(alphaVec(lp1),deltaVec(lp2));
    end
end

%Plot
surf(X,Y,Z,abs(fVals));
shading interp;
