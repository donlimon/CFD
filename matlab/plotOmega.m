clear all;
delimiterIn = ' ';
% DATA FILES
x_file = '../data/grid-x.dat';
y_file = '../data/grid-y.dat';
omega_file = '../data/data-omega.dat';
% IMPORT FILES
x = importdata(x_file,delimiterIn);
y = importdata(y_file,delimiterIn);
omega = importdata(omega_file,delimiterIn);
% PLOT
contourf(x,y,omega);
colorbar;
