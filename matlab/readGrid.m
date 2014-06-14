function [x y] = readGrid(dataDir,ngp)
	xFileHandle = fopen([dataDir 'grid-x.bin']);
	yFileHandle = fopen([dataDir 'grid-y.bin']);
	
	x = fread(xFileHandle,[ngp ngp],'double');
	y = fread(yFileHandle,[ngp ngp],'double');
	
	fclose(xFileHandle);
	fclose(yFileHandle);
end
