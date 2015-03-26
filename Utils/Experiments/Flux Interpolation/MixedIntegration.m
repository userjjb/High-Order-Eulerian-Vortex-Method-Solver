clear all

M=8;
N=11;
[Qx,Qw]= GLquad(M);
[nd,w] = gauss(N);

mm=elim(Qx(1:M)',Qx(1:M)',[1 3 2]);
nn=elim(nd(1:N)',nd(1:N)',[1 3 2]);
Lag2= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,mm(nv,:,:)),bsxfun(@minus,Qx(nv),mm(nv,:,:))),3);
Lag= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,nn(nv,:,:)),bsxfun(@minus,nd(nv),nn(nv,:,:))),3);
dift= @(x,n) sum(1./bsxfun(@minus,x,nd([1:n-1,n+1:N])));
dLag= @(x,n) Lag(x,n).*dift(x,n);

s = @(x) sin(pi*x);
    c2 = 60/100;
    b=-.2;
g = @(x) exp(-(x-b).^2/(2*c2));
s_g=@(x) s(x).*g(x);
    
interp_s=@(x) s(nd')*Lag(x,1:N);
interp_g=@(x) g(Qx')*Lag2(x,1:M);
interp_s_g=@(x) s_g(nd')*Lag(x,1:N);
interp_sg=@(x) interp_s(x).*interp_g(x);

n=2;
exact=@(x) s_g(x).*dLag(x-eps(x),n);
stiff_s_g=@(x) interp_s_g(x).*dLag(x-eps(x),n);
stiff_sg=@(x) interp_sg(x).*dLag(x-eps(x),n);
    
E=integral(exact,-1,1,'RelTol',1e-14,'AbsTol',1e-16)
I=integral(stiff_sg,-1,1,'RelTol',1e-14,'AbsTol',1e-16)
G=stiff_s_g(nd')*w'

for i=1:M
    for j=1:N
            trip=@(x) Lag2(x,i).*Lag(x,j).*dLag(x,n);
            W(i,j)=integral(trip,-1,1,'RelTol',1e-14,'AbsTol',1e-16);
    end
end
J=interp_g(Qx')*(W*interp_s(nd')')

xx=-1:.01:1;
plot(xx,s_g(xx),'k')
hold on
plot(xx,interp_s_g(xx),'r')
plot(xx,interp_sg(xx),'g')