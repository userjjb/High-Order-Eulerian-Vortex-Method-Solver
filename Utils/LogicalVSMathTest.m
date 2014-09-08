s = 4000;
tic
not(toeplitz(mod(0:s,2),0:s<0));
toc

tic
ones(s+1)-(toeplitz(mod(0:s,2),0:s<0));
toc