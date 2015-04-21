clear all
N=8;
del=0.01;
r=@(x) x./(x.^2+del^2).^(3/2);
A=0;
B=1;

[Qx, Qw] = GenOrthog(N,r,A,B,1);

%In case you're tempted to think just sticking a x=0 abcissa in the front
%makes things better:
% [Qx3, Qw3]= GLquad(N-1);
% Qx3=[0; (Qx3+1)/2];
% nn3=elim(Qx3(1:N)',Qx3(1:N)',[1 3 2]);
% Lag3= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,nn3(nv,:,:)),bsxfun(@minus,Qx3(nv),nn3(nv,:,:))),3);
% for node=1:N
%     %We'll likely throw some convergence warnings here
%     Qw3(1,node)=integral(@(x) Lag3(x,node),A,B,'RelTol',1e-14,'AbsTol',1e-16);
% end

[Qx2, Qw2]= GLquad(N);
Qx2=(Qx2+1)/2;

%fun = @(x) sin(pi*x-pi/2)+1;
%fun = @(x) max(0, exp(-(x-.5).^2/.2));
fun = @(x) 1-(pi^2*x.^2)/2+(pi^4*x.^4)/24-(pi^6*x.^6)/720; 

nn=elim(Qx2',Qx2',[1 3 2]);
Lag= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,nn(nv,:,:)),bsxfun(@minus,Qx2(nv),nn(nv,:,:))),3);
interp_f=@(x) fun(Qx2')*Lag(x,1:N);

d=0;

E=integral(@(x) fun(x).*r(x+d),0,1,'RelTol',1e-14,'AbsTol',1e-16)
I=integral(@(x) interp_f(x).*r(x+d),0,1,'RelTol',1e-14,'AbsTol',1e-16)
Q=(fun(Qx2').*r(Qx2'+d))*Qw2'/2
QJ=(fun(Qx').*r(Qx'+d))*Qw'
%QJi=(interp_f(Qx').*r(Qx'+d))*Qw'

%For del 0.01
Qx3=[0.0102416656563158,0.0717853629985321,0.270109383343779,0.367444573443785,0.547385302479117,0.698468748230208,0.833622761210289,0.970714987883877];
Qw3=[0.0198977035216784,0.130830119885271,0.246726289811717,0.0140463368447291,0.257974748739917,0.0658727025304273,0.186175240124866,0.0784768585413942];
Q3=fun(Qx3).*r(Qx3)*Qw3'
%Q3=(fun(Qx3').*r(Qx3'))*Qw3'