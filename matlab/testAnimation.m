clear all;

readSettings;
[vor1 vor2] = readVortexData(tsmax);
animateVortexCenter(vor1,vor2,dx,tsmax);
