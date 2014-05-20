clear all;
% Get data from binary files
xFile = fopen('../data/grid-x.bin');
yFile = fopen('../data/grid-y.bin');
uFile = fopen('../data/data-u.bin');
vFile = fopen('../data/data-v.bin');
mFile = fopen('../data/data-vmag.bin');
pFile = fopen('../data/data-p.bin');

xFileInfo = dir('../data/grid-x.bin');
yFileInfo = dir('../data/grid-y.bin');
uFileInfo = dir('../data/data-u.bin');
vFileInfo = dir('../data/data-v.bin');
mFileInfo = dir('../data/data-vmag.bin');
pFileInfo = dir('../data/data-p.bin');

x = fread(xFile,[sqrt(xFileInfo.bytes/8) sqrt(xFileInfo.bytes/8)],'double');
y = fread(yFile,[sqrt(yFileInfo.bytes/8) sqrt(yFileInfo.bytes/8)],'double');
u = fread(uFile,[sqrt(uFileInfo.bytes/8) sqrt(uFileInfo.bytes/8)],'double');
v = fread(vFile,[sqrt(vFileInfo.bytes/8) sqrt(vFileInfo.bytes/8)],'double');
m = fread(mFile,[sqrt(mFileInfo.bytes/8) sqrt(mFileInfo.bytes/8)],'double');
p = fread(pFile,[sqrt(pFileInfo.bytes/8) sqrt(pFileInfo.bytes/8)],'double');

% Plot data
figure;
vel_magnitude = contourf(x,y,m);
colorbar;
figure;
vel_vector = quiver(x,y,u,v);
figure;
pressure = contourf(x,y,p);
colorbar;
