clear all
xx= -1:.00001:1;
q= .5;
core= 0.1;
yy= (cos(pi*xx)/2+.5).*(q-xx)./abs(q-xx).^3;
yyRM= (cos(pi*xx)/2+.5).*(q-xx)./((q-xx).^2+core^2).^(3/2);
hold on
axis([-1 1 min(yyRM) max(yyRM)])
plot(xx,yy,'b')
plot(xx,yyRM,'r')
text(.5,1,num2str(find(max(yyRM)==yyRM)))%Resolved discontinuity location (displaced from actual singularity, depending on core radius)