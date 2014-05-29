clear all;
dataPath = '../data/';
fileInfo = dir([dataPath 'grid-x.bin']);
NGP = sqrt(fileInfo.bytes/8);

% Get data from binary files
xFile = fopen([dataPath 'grid-x.bin']);
yFile = fopen([dataPath 'grid-y.bin']);

x = fread(xFile,[NGP NGP],'double');
y = fread(yFile,[NGP NGP],'double');

uFiles=dir([dataPath 'data-u-*.bin']);
for i=1:length(uFiles)
	uFile(i) = fopen([dataPath uFiles(i).name]);
    u(i) = fread(uFile(i),[NGP NGP],'double');
    fclose(uFile(i));
end

vFiles=dir([dataPath 'data-v-*.bin']);
for i=1:length(vFiles)
	vFile(i) = fopen([dataPath vFiles(i).name]);
    v(i) = fread(vFile(i),[NGP NGP],'double');
    fclose(vFile(i));
end

mFiles=dir([dataPath 'data-vmag-*.bin']);
for i=1:length(mFiles)
	mFile(i) = fopen([dataPath mFiles(i).name]);
    m(i) = fread(mFile(i),[NGP NGP],'double');
    fclose(mFile(i));
end

pFiles=dir([dataPath 'data-p-*.bin']);
for i=1:length(pFiles)
	pFile(i) = fopen([dataPath pFiles(i).name]);
    p(i) = fread(pFile(i),[NGP NGP],'double');
    fclose(pFile(i));
end

% Plot data
%figure;
%vel_magnitude = contourf(x,y,m);
%colorbar;
%figure;
%vel_vector = quiver(x,y,u,v);
%figure;
%pressure = contourf(x,y,p);
%colorbar;
