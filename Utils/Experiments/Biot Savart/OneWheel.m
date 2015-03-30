a=0.1;
b=0.1;
g=@(x,y) 15*exp(-((x-.5).^2/a+(y-.5).^2/b));
r=@(x,y,i,j) 1./((i+x).^2+(j+y).^2);

N=10;
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
for i=1:20
    for j=1:20
        duo=@(x,y) g(x,y).*r(x,y,i,j);
        iduo=@(x,y) interp_g(x,y).*r(x,y,i,j);
        E(i,j)=integral2(duo,0,1,0,1);
        I(i,j)=integral2(iduo,0,1,0,1);
        for k=1:N
            for m=1:N
                R(k,m)=r(nd(k),nd(m),i,j);
            end
        end
        Q(i,j)=Qw/4*(G.*R)*Qw';
        OW(i,j)=q/((i+.5)^2+(j+.5)^2);
    end
end