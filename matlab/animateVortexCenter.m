function [] = animateVortexCenter(vor1,vor2,dx,tsmax)
	figure;
	handle = plot(vor1(1).i*dx,vor1(1).j*dx,'kx');
	axis([0 1 0 1]);
	hold on;
	plot(vor2(1).i*dx,vor2(1).j*dx,'kx');
	for i=2:tsmax
		plot([vor1(i).i*dx vor1(i-1).i*dx],[vor1(i).j*dx vor1(i-1).j*dx],'-k');
		plot([vor2(i).i*dx vor2(i-1).i*dx],[vor2(i).j*dx vor2(i-1).j*dx],'-k');
		%pause(0.001);
	end
	hold off;
end
