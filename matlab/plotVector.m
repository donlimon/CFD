clear all;
delimiterIn = ' ';
% DATA FILES
x_file = '../data/grid-x.dat';
y_file = '../data/grid-y.dat';
u_file = '../data/data-u.dat';
v_file = '../data/data-v.dat';
% IMPORT FILES
x = importdata(x_file,delimiterIn);
y = importdata(y_file,delimiterIn);
u = importdata(u_file,delimiterIn);
v = importdata(v_file,delimiterIn);
% PLOT
quiver(x,y,u,v);
