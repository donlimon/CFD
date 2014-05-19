clear all;
delimiterIn = ' ';
% DATA FILES
x_file = '../data/grid-x.dat';
y_file = '../data/grid-y.dat';
phi_file = '../data/data-phi.dat';
% IMPORT FILES
x = importdata(x_file,delimiterIn);
y = importdata(y_file,delimiterIn);
phi = importdata(phi_file,delimiterIn);
% PLOT
contourf(x,y,phi);
colorbar;
