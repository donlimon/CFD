clear all;
dataDir = '../data/';
readSettings;
[x y] = readGrid(dataDir,ngp);
omega = readVorticity(dataDir,ngp,2500,tsmax);

handle = figure();
contourf(x,y,omega);
[vor1 vor2] = readVortexData(dataDir,tsmax);
plotVortexCenter(handle,vor1,vor2,2500,dx);
