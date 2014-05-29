clear all;
readSettings;
[x y] = readGrid(ngp);

handle = figure();
plotVelmag(handle,x,y,ngp,500,tsmax);
[vor1 vor2] = readVortexData(tsmax);
plotVortexCenter(handle,vor1,vor2,500,dx);
