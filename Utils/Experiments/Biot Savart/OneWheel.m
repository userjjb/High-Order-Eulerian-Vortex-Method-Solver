clear all

a=0.3;
b=0.7;
dx=.1; 
dy=.2;
g=@(x,y) 15*exp(-((x-dx).^2/a+(y-dy).^2/b));%Center is at (dx,dy)
r=@(x,y,i,j) 1./((i+x).^2+(j+y).^2);

N=7;
[Qx,Qw] = GLquad(N);

nd=(Qx+1)/2;
nn=elim(nd(1:N)',nd(1:N)',[1 3 2]);
Lag= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,nn(nv,:,:)),bsxfun(@minus,nd(nv),nn(nv,:,:))),3);
for i=1:N
    for j=1:N
        G(i,j)=g(nd(i),nd(j));
    end
end
interp_g=@(x,y) reshape(diag(Lag(reshape(x,1,[]),1:N)'*(G*Lag(reshape(y,1,[]),1:N)))',size(x));

q=integral2(g,0,1,0,1);
qq=integral2(interp_g,0,1,0,1);
for i=0:20
    ii=i+1;
    for j=0:20
        jj=j+1;
        duo=@(x,y) g(x,y).*r(x,y,i,j);
        iduo=@(x,y) interp_g(x,y).*r(x,y,i,j);
        E(ii,jj)=integral2(duo,0,1,0,1);
        I(ii,jj)=integral2(iduo,0,1,0,1);
        for k=1:N
            for m=1:N
                R(k,m)=r(nd(k),nd(m),i,j);
            end
        end
        Q(ii,jj)=Qw/4*(G.*R)*Qw';
        OW(ii,jj)=qq/((i+.332)^2+(j+.429)^2);
    end
end

difEQ=E-Q;
difEI=E-I;
difIQ=I-Q;