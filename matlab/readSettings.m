setFile = fopen('../data/settings.bin','r');

simtype = char(fread(setFile,1,'char*1'));
ngp = fread(setFile,1,'int');
tsmax = fread(setFile,1,'int');
lref = fread(setFile,1,'double');
dt = fread(setFile,1,'double');
dx = fread(setFile,1,'double');

if simtype=='t'
	Re = fread(setFile,1,'double');
	uref = fread(setFile,1,'double');
elseif simtype=='v'
	a = fread(setFile,1,'double');
	lratio = fread(setFile,1,'double');
	Re = fread(setFile,1,'double');
	nu = fread(setFile,1,'double');
end

fclose(setFile);
clear setFile;
