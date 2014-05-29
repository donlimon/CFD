function [] = plotVelmag(handle,x,y,ngp,ts,tsmax)
	digits = 1+floor(log10(abs(tsmax)));
	figure(handle);
	mFile = fopen(['../data/data-vmag-' sprintf(['%0' num2str(digits) 'd'],ts) '.bin']);
	m = fread(mFile,[ngp ngp],'double');
	hold on;
	contourf(x,y,m);
	colorbar;
	hold off;
end
