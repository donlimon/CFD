function [x y] = readGrid(ngp)
	xFile = fopen('../data/grid-x.bin');
	yFile = fopen('../data/grid-y.bin');
	
	x = fread(xFile,[ngp ngp],'double');
	y = fread(yFile,[ngp ngp],'double');
end
