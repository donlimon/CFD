function [vortex1 vortex2] = readVortexData(tsmax)
	% Get data from binary files
	vorFile1 = fopen('../data/data-vortex-1.bin');
	vorFile2 = fopen('../data/data-vortex-2.bin');

	%allocate memory
	vortex1(tsmax) = struct('i',[], 'j',[], 'ts',[]);
	vortex2(tsmax) = struct('i',[], 'j',[], 'ts',[]);

	%read vortex data
	for i=1:tsmax
		%vortex #1
		vortex1(i).i = fread(vorFile1,1,'int');
		vortex1(i).j = fread(vorFile1,1,'int');
		vortex1(i).ts = fread(vorFile1,1,'int');
		%vortex #2
		vortex2(i).i = fread(vorFile2,1,'int');
		vortex2(i).j = fread(vorFile2,1,'int');
		vortex2(i).ts = fread(vorFile2,1,'int');
	end

	fclose(vorFile1);
	fclose(vorFile2);
end
