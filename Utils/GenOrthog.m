clear all

N=4;

del=0.01;
w=@(x) x./(x.^2+del^2).^(3/2);

p{1}=@(x) 1;
c2= integral(@(x) x.*p{1}(x).^2.*w(x),0,1,'RelTol',1e-14,'AbsTol',1e-16)/integral(@(x) p{1}(x).^2.*w(x),0,1,'RelTol',1e-14,'AbsTol',1e-16);
p{2}=@(x) (x-c2).*p{1}(x);

I2(1)= integral(@(x) p{1}(x).^2.*w(x),0,1);
for iter=3:N+1
    I2(iter-1)= integral(@(x) p{iter-1}(x).^2.*w(x),0,1,'RelTol',1e-14,'AbsTol',1e-16);
    c2= integral(@(x) x.*p{iter-1}(x).^2.*w(x),0,1,'RelTol',1e-14,'AbsTol',1e-16)/I2(iter-1);
    c1= I2(iter-1)/I2(iter-2);
    p{iter}=@(x) (x-c2).*p{iter-1}(x) - c1*p{iter-2}(x);
end

orthRoot(2,1)= fzero(p{2},[0 1]);
for poly=3:N+1
    for numroot=1:poly-1
        temp= [0, orthRoot(poly-1,:), 1];
        orthRoot(poly,numroot)= fzero(p{poly},[temp(numroot), temp(numroot+1)]);
    end
end

Qx=orthRoot(N,:)';
nn=elim(Qx(1:N)',Qx(1:N)',[1 3 2]);
Lag= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,nn(nv,:,:)),bsxfun(@minus,Qx(nv),nn(nv,:,:))),3);

for node=1:N
    Qw(1,node)=integral(@(x) Lag(x,node),0,1,'RelTol',1e-14,'AbsTol',1e-16);
end

fun = @(x) sin(pi*x-pi/2)+1;
del=0.01;
r=@(x) x./(x.^2+del^2).^(3/2);
QJ=(fun(Qx').*r(Qx'))*Qw'
E=integral(@(x) fun(x).*r(x),0,1,'RelTol',1e-14,'AbsTol',1e-16)

[Qx2, Qw2]= GLquad(N);
Qx2=(Qx2+1)/2;
Q=(fun(Qx2').*r(Qx2'))*Qw2'/2