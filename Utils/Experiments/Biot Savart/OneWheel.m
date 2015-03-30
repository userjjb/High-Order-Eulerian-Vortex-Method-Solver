a=0.3;
b=0.7;
dx=.1; 
dy=.2;
g=@(x,y) 15*exp(-((x-dx).^2/a+(y-dy).^2/b));%Center is at (dx,dy)
r=@(x,y,i,j) 1./((i+x).^2+(j+y).^2);

N=8;
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
    for j=1:20
        duo=@(x,y) g(x,y).*r(x,y,i,j);
        iduo=@(x,y) interp_g(x,y).*r(x,y,i,j);
        E(ii,j)=integral2(duo,0,1,0,1);
        I(ii,j)=integral2(iduo,0,1,0,1);
        for k=1:N
            for m=1:N
                R(k,m)=r(nd(k),nd(m),i,j);
            end
        end
        Q(ii,j)=Qw/4*(G.*R)*Qw';
        OW(ii,j)=qq/((i+.332)^2+(j+.429)^2);
    end
end

difEOW=E-OW;
difEQ=E-Q;
difQOW=Q-OW;
rrr=q./E;
rrrs=sqrt(rrr);
[Y2, X2]= meshgrid(1:20,0:20);
rr=(X2+.332).^2+(Y2+.429).^2;
rrd=rr-rrr;
[X, Y]= meshgrid(0:.01:1,0:.01:1);
%surf(X,Y,g(X,Y))
norm(rrd)