readSettings;
dataDir='../data/';
[x y] = readGrid(dataDir,ngp);
figure();
for i=50:50:5000
	p = readPressure(dataDir,ngp,i,tsmax);
	contourf(x,y,p);
	pause(0.1);
end
