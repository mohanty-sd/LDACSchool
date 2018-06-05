%% Test SKYPLOT
alphaVec = 0:0.05:(2*pi);
deltaVec = 0:0.05:pi;

%Test function
fHandle = @(x,y) sin(2*x)*cos(2*y);

%Run skyplot
skyplot(alphaVec,deltaVec,fHandle);
axis equal;