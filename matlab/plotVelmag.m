clear all;
delimiterIn = ' ';
% DATA FILES
x_file = '../data/grid-x.dat';
y_file = '../data/grid-y.dat';
vel_file = '../data/data-velmag.dat';
% IMPORT FILES
x = importdata(x_file,delimiterIn);
y = importdata(y_file,delimiterIn);
vel = importdata(vel_file,delimiterIn);
% PLOT
contourf(x,y,vel);
colorbar;
