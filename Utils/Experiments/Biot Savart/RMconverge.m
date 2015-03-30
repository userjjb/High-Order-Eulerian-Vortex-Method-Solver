clear all

K=@(x,d) x./(x.^2 +d^2).^(3/2);
KK=@(x) 1./x.^2;
    c2 = 5/100;
w=@(x) exp(-(x+.1).^2/(2*c2));

duo=@(x) K(x+1,.1).*w(x);
q=integral(duo,-1,1,'RelTol',1e-14,'AbsTol',1e-16)

N=9;
[Qx, Qw]= GLquad(N);
pp=elim(Qx(1:N)',Qx(1:N)',[1 3 2]);
Lag3= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,pp(nv,:,:)),bsxfun(@minus,Qx(nv),pp(nv,:,:))),3);

g=duo(Qx')*Qw'
interp=@(x) duo(Qx')*Lag3(x,1:N);

xx=-1:.01:1;
plot(xx,duo(xx),'k')
hold on
plot(xx,interp(xx),'r')