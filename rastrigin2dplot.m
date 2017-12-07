%% 2D Rastrigin Function
% Plots the 2D Rastrigin function, which is one of the standard benchmark
% fitness functions used for characterizing the performance of global
% minimization function.

%%
% Authors:
% Soumya D. Mohanty, Dec 2017
% First code
%
%

%%
% Grid of points along each axis
x = -5:.1:5;
y = -5:.1:5;
%%
% Convert to 2D array of grid points 
% X: each row is x
% Y: each column is y
[X,Y]=meshgrid(x,y);
%%
% f(X)+g(Y) = Z where Z(i,j) = f(X(i,j))+g(Y(i,j)) = f(x(j))+g(y(i))
figure;
surf(X,Y,X.^2+Y.^2-10*cos(2*pi*X)-10*cos(2*pi*Y)+2*10)
xlabel('x');
ylabel('y');
title('Rastrigin function of 2 variables: \Sigma_i (x_i^2-10*cos(2\pi x_i)+10)');
shading interp