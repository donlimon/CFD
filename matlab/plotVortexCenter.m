function [] = plotVortexCenter(handle,vortex1,vortex2,ts,dx)
	figure(handle);
	hold on;
	plot(vortex1(ts).i*dx,vortex1(ts).j*dx,'kx','MarkerSize', 30);
	plot(vortex2(ts).i*dx,vortex2(ts).j*dx,'kx','MarkerSize', 30);
	hold off;
end
