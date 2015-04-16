clear all

N=5;
[Qx,Qw] = GLquad(N);
nd=(Qx+1)/2;

a=0.3;
b=0.3;
dx=0.8; 
dy=0.5;
w=@(x,y) 15*exp(-((x-dx).^2/a+(y-dy).^2/b));%Center is at (dx,dy)

Tx=.3;
Ty=nd(N)-1;%
K=@(x,y,i,j) ((x+i)-Tx)./((i+x-Tx).^2+(j+y-Ty).^2).^(3/2);
%----------
mid=0;
r=@(x) 1./((x-Tx).^2+(mid-Ty).^2).^(3/2);
[QxSeed,QwSeed] = GenOrthog(N,r,0,1);
for Sn=1:N
    r=@(x) (QxSeed(Sn)-Tx)./((x-Tx).^2+(QxSeed(Sn)-Ty).^2).^(3/2);
    [QxJ, QwJ] = GenOrthog(N,r,0,1);
    QxJM(:,Sn)=QxJ;
    QwJM(Sn,:)=QwJ;
end

for Sn=1:N
    r=@(x) 1./((nd(Sn)-Tx).^2+(x-Ty).^2).^(3/2);
    [QxJ, QwJ] = GenOrthog(N,r,0,1);
    QxJMe(:,Sn)=QxJ;
    QwJMe(Sn,:)=QwJ;
end

%----------
mid2=0;
r=@(x) 1./((x-Tx).^2+(mid2-Ty).^2).^(3/2);
[QxJ,QwJ] = GenOrthog(N,r,0,1);
%----------
nn=elim(nd(1:N)',nd(1:N)',[1 3 2]);
Lag= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,nn(nv,:,:)),bsxfun(@minus,nd(nv),nn(nv,:,:))),3);
for i=1:N
    for j=1:N
        W(i,j)=w(nd(j),nd(i));
    end
end
interp_w=@(x,y) reshape(diag(Lag(reshape(x,1,[]),1:N)'*(W'*Lag(reshape(y,1,[]),1:N)))',size(x));

for i=0:20
    ii=i+1;
    for j=0:20
        jj=j+1;
        duo=@(x,y) w(x,y).*K(x,y,i,j);
        iduo=@(x,y) interp_w(x,y).*K(x,y,i,j);
        E(jj,ii)=integral2(duo,0,1,0,1);
        I(jj,ii)=integral2(iduo,0,1,0,1);
        for k=1:N
            for m=1:N
                R(m,k)=K(nd(k),nd(m),i,j);
            end
        end
        Q(jj,ii)=Qw/4*(W.*R)*Qw';
    end
end
%-----
for i=1:N
    for j=1:N
        WJi(i,j)=interp_w(QxJM(i,j),QxSeed(j));
    end
end

for k=1:N
    for m=1:N
        R(k,m)=K(QxJM(k,m),QxSeed(m),0,0);
    end
end
QJ=sum(QwJM'.*(WJi.*R))*QwSeed';
E(1,1)-QJ
%-----
for i=1:N
    for j=1:N
        WTJi(j,i)=interp_w(QxJ(i),QxJ(j));
    end
end

for k=1:N
    for m=1:N
        RT(m,k)=K(QxJ(k),QxJ(m),0,0);
    end
end
QTJ=QwJ*(WTJi.*RT)*QwJ';
E(1,1)-QTJ
%-----
for i=1:N
    for j=1:N
        WeJi(i,j)=interp_w(nd(j),QxJMe(i,j));
    end
end

for k=1:N
    for m=1:N
        Re(k,m)=K(nd(m),QxJMe(k,m),0,0);
    end
end
QJe=sum(QwJMe'.*(WeJi.*Re))*Qw'/2;
E(1,1)-QJe
%-----
for i=1:N
    for j=1:N
        WemJi(i,j)=interp_w(nd(j),QxSeed(i));
    end
end

for k=1:N
    for m=1:N
        Rem(k,m)=K(nd(m),QxSeed(k),0,0);
    end
end
QJem=QwSeed*(WemJi.*Rem)*Qw'/2;
E(1,1)-QJem

difEQ=E-Q;
difEI=E-I;
difIQ=I-Q;