clear all

M=9;
[Qx,Qw]= GLquad(M);

nn=elim(Qx(1:M)',Qx(1:M)',[1 3 2]);
Lag2= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,nn(nv,:,:)),bsxfun(@minus,Qx(nv),nn(nv,:,:))),3);

s = @(x) sin(pi*x);
    c2 = 20/100;
    b=-.2;
g = @(x) exp(-(x-b).^2/(2*c2));
s_g=@(x) s(x).*g(x);

for N=3:20
    [nd,w] = gauss(N);
    
    nn=elim(nd(1:N)',nd(1:N)',[1 3 2]);
    Lag= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,nn(nv,:,:)),bsxfun(@minus,nd(nv),nn(nv,:,:))),3);
    dift= @(x,n) sum(1./bsxfun(@minus,x,nd([1:n-1,n+1:N])));
    dLag= @(x,n) Lag(x,n).*dift(x,n);
    
    interp_s=@(x) s(nd')*Lag(x,1:N);
    interp_g=@(x) g(nd')*Lag(x,1:N);
    interp_s_g=@(x) s_g(nd')*Lag(x,1:N);
    interp_sg=@(x) interp_s(x).*interp_g(x);
    
    n=floor(N/2);
    exact=@(x) s_g(x).*dLag(x-eps(x),n);
    stiff_s_g=@(x) interp_s_g(x).*dLag(x-eps(x),n);
    stiff_sg=@(x) interp_s(x).*interp_g(x).*dLag(x-eps(x),n);
    
    Q(N,3)=integral(exact,-1,1,'RelTol',1e-14,'AbsTol',1e-16);
    Q(N,4)=integral(stiff_s_g,-1,1,'RelTol',1e-14,'AbsTol',1e-16);
    Q(N,5)=integral(stiff_sg,-1,1,'RelTol',1e-14,'AbsTol',1e-16);
    
    Q(N,6)=stiff_s_g(Qx')*Qw';
    Q(N,7)=stiff_sg(Qx')*Qw';
        
    for i=1:M
        duo=@(x) Lag2(x,i).*dLag(x,n);
        W(i)=integral(duo,-1,1,'RelTol',1e-14,'AbsTol',1e-15);
    end
    Q(N,8)=interp_s_g(Qx')*W';
    
    for i=1:M
        for j=1:M
                trip=@(x) Lag2(x,i).*Lag2(x,j).*dLag(x,n);
                W2(i,j)=integral(trip,-1,1,'RelTol',1e-14,'AbsTol',1e-15);
        end
    end
    Q(N,9)=(interp_s(Qx')*W2)*interp_g(Qx')';
    Q(N,13)=(N-2);
end
Q(:,1)=Q(:,3)-Q(:,4);
Q(:,2)=Q(:,3)-Q(:,5);

Q(:,10)=Q(:,4)-Q(:,6);
Q(:,11)=Q(:,4)-Q(:,8);
Q(:,12)=Q(:,5)-Q(:,9);