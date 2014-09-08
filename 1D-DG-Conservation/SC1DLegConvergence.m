h = waitbar(0,'Initializing waitbar...');
prog=0;
for N=1:9
    for KK=1:8
        K=2^KK;
        prog = prog+(N*K);
        waitbar(prog/22950,h,'Calculating')
        saved(N,KK) = SC1DLeg(N,K);
    end
end
close(h)