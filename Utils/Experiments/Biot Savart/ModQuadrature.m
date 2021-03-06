clear all

N=4;
[Qx, Qw]= GLquad(N);

%Local vorticity functions
%fun = @(x) 1-x.^2/2+x.^4/24-x.^6/720+x.^8/40320-x.^10/3628800; %Taylor trig
%fun = @(x) 1-(pi^2*x.^2)/2+(pi^4*x.^4)/24-(pi^6*x.^6)/720;     %Taylor trig
%fun = @(x) max(0, exp(-(x-.5).^2/.04));                              %Gaussian curve
%fun = @(x) sin(pi*x-pi/2)+1;       %Smooth trig
%fun = @(x) max(0,sin(pi*(x-.1)));       %Cleaved trig
%fun = @(x) heaviside(x+0).*exp(-1./(1-(2*(x-.5)).^2)).*heaviside(1-x); %Mollifier
del=0.01;
r=@(x) x./(x.^2+del^2).^(3/2); %RM
%r=@(x) ((x.^2+2.5*del^2).*x)./((x+d).^2+del^2).^(5/2); %WL

nd=(Qx+1)/2;
nn=elim(nd(1:N)',nd(1:N)',[1 3 2]);
Lag= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,nn(nv,:,:)),bsxfun(@minus,nd(nv),nn(nv,:,:))),3);
interp_fun=@(x) fun(nd')*Lag(x,1:N);

M=10;
[Qx2, Qw2]= GLquad(M);

duo=@(x) fun(x).*r(x);
iduo=@(x) interp_fun(x).*r(x);
E=integral(duo,0,1,'RelTol',1e-14,'AbsTol',1e-16)
I=integral(iduo,0,1,'RelTol',1e-14,'AbsTol',1e-16)
Q=(fun(nd').*r(nd'))*Qw'/2
%Q2=(fun(nd(1:end-1)')./nd(1:end-1)'.^2)*Qw(1:end-1)'/2
E-I
I-Q
xx=0:.001:1;