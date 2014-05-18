clear all;
delimiterIn = ' ';
% DATA FILES
x_file = '../data/grid-x.dat';
y_file = '../data/grid-y.dat';
p_file = '../data/data-p.dat';
% IMPORT FILES
x = importdata(x_file,delimiterIn);
y = importdata(y_file,delimiterIn);
p = importdata(p_file,delimiterIn);
% PLOT
contourf(x,y,p);
colorbar;
