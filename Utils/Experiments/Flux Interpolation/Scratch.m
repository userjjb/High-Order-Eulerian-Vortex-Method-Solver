clear all

P=10;
[Qx,Qw]= GLquad(P);
pp=elim(Qx(1:P)',Qx(1:P)',[1 3 2]);
Lag3= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,pp(nv,:,:)),bsxfun(@minus,Qx(nv),pp(nv,:,:))),3);

R=10;
[Qx2,Qw2]= LGLquad(R);
rr=elim(Qx2(1:R)',Qx2(1:R)',[1 3 2]);
Lag4= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,rr(nv,:,:)),bsxfun(@minus,Qx2(nv),rr(nv,:,:))),3);

s = @(x) sin(pi*x);
    c2 = 40/100;
    b=-.2;
g = @(x) exp(-(x-b).^2/(2*c2));
s_g=@(x) s(x).*g(x);

for N=3:15
    [nd,w] = GLquad(N);
    M=N;
    [nd2,w2] = GLquad(M);
    
    %Fully vectorized for both x and nv
    nn=elim(nd(1:N)',nd(1:N)',[1 3 2]);
    Lag= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,nn(nv,:,:)),bsxfun(@minus,nd(nv),nn(nv,:,:))),3);
    nn2=elim(nd2(1:M)',nd2(1:M)',[1 3 2]);
    Lag2= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,nn2(nv,:,:)),bsxfun(@minus,nd2(nv),nn2(nv,:,:))),3);
    dLag= @(x,nv) Lag(x,nv).*sum(1./bsxfun(@minus,x,nn(nv,:,:)),3);
    
    interp_s=@(x) s(nd')*Lag(x,1:N);
    interp_g=@(x) g(nd2')*Lag2(x,1:M);
    interp_s_g=@(x) s_g(nd2')*Lag2(x,1:M);
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
        
    for i=1:P
        duo=@(x) Lag3(x,i).*dLag(x,n);
        W(i)=integral(duo,-1,1,'RelTol',1e-14,'AbsTol',1e-15);
    end
    Q(N,8)=interp_s_g(Qx')*W';
    
    for i=1:P
        for j=1:R
                trip=@(x) Lag3(x,i).*Lag4(x,j).*dLag(x,n);
                W2(i,j)=integral(trip,-1,1,'RelTol',1e-14,'AbsTol',1e-15);
        end
    end
    Q(N,9)=(interp_s(Qx')*W2)*interp_g(Qx2')';
    
    Q(N,15)=(N-1)+(M-1);
end
Q(:,1)=Q(:,3)-Q(:,4);
Q(:,2)=Q(:,3)-Q(:,5);

Q(:,11)=Q(:,4)-Q(:,6);
Q(:,12)=Q(:,4)-Q(:,8);
Q(:,13)=Q(:,5)-Q(:,9);
Q(:,14)=Q(:,3)-Q(:,9);