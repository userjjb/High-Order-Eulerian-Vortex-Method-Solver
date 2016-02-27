ccA=rand(6,1,21600);
ccB=rand(1,6,21600);
ccC=permute(ccA,[1 3 2]); %[1 3 2]
ccD=permute(ccB,[2 3 1]); %[2 3 1]

for it=1:10000
    ccE= mtimesx(ccB,ccA);
    ccF= sum(ccC.*ccD);
end