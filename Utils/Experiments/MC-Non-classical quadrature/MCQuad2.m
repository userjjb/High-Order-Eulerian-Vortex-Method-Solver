clear all
syms xx yy

del=0.01;
N=5;

[Qx,Qw] =GLquad(N);
nd=(Qx+1)/2;
[ndx, ndy] = meshgrid(nd,nd);

[Qx2,Qw2] =GLquad(N+1);
nd2=(Qx2+1)/2;
[ndx2, ndy2] = meshgrid(nd2,nd2);

Ty=0;
for nx=1:N
Tx=nd(nx);

%Test cases
k=@(x,y) (y-Ty)./((x-Tx).^2+(y-Ty).^2+del^2).^(3/2);

a=0.2;
b=0.3;
dx=0.6; 
dy=0.6;
w=@(x,y) exp(-((x-dx).^2/a+(y-dy).^2/b));
T1{nx}=@(x,y) w(x,y).*k(x,y);
E(1,nx)=integral2(T1{nx},0,1,0,1);
Q(1,nx)=Qw*T1{nx}(ndx,ndy)*Qw'/4;

a=0.05;
b=0.05;
dx=0.5; 
dy=0.5;
w=@(x,y) exp(-((x-dx).^2/a+(y-dy).^2/b));
T2{nx}=@(x,y) w(x,y).*k(x,y);
E(2,nx)=integral2(T2{nx},0,1,0,1);
Q(2,nx)=Qw*T2{nx}(ndx,ndy)*Qw'/4;

a=0.06;
b=0.15;
dx=0.1; 
dy=0.6;
w=@(x,y) exp(-((x-dx).^2/a+(y-dy).^2/b));
T3{nx}=@(x,y) w(x,y).*k(x,y);
E(3,nx)=integral2(T3{nx},0,1,0,1);
Q(3,nx)=Qw*T3{nx}(ndx,ndy)*Qw'/4;
end
err=abs(E-Q);

%Begin MC method
p=nchoosek(N+2,N);

xe=0;
ye=0;
for n=1:N
    xe=[xe, 0:n];
    ye=[ye, n:-1:0];
end

%Old slower symbolic math way
% inc=1;
% for n=0:N
%     for i=0:n
%         R(inc)= xx^i*yy^(n-i);
%         inc=inc+1;
%     end
% end
% 
% dM=det(M);
% for B=1:p
%     Mb=[M(1:B-1,:);R;M(B+1:end,:)];
%     dB=det(Mb);
%     QwC(B)=integral2(matlabFunction(dB/dM),0,1,0,1);
% end

QuadDet=bsxfun(@power,reshape(ndx2,[],1),xe).*bsxfun(@power,reshape(ndy2,[],1),ye);

inc=1;
inc2=1;
C=1E6;
skip=0;
best=0;
while inc<1000   
    while C>4E5
        Qx=rand(p,2);
        M=bsxfun(@power,Qx(:,1),xe).*bsxfun(@power,Qx(:,2),ye);
        C=cond(M);
    end
    conder(inc2)=C;
    C=1E6;
    for B=1:p
        for nodes=1:(N+1)^2
            dB(nodes)=det([M(1:B-1,:);QuadDet(nodes,:);M(B+1:end,:)]);
        end
        Qw(B)=Qw2*reshape(dB,N+1,[])*Qw2'/4;
%         if Qw(B)<0
%             skip=1;
%             break
%         end
    end
%     if skip
%         skip=0;
%         continue
%     end
    Qw=Qw/det(M);
    
    for nx=1:1
        errB(1,nx)=T1{nx}(Qx(:,1),Qx(:,2))'*Qw';
        errB(2,nx)=T2{nx}(Qx(:,1),Qx(:,2))'*Qw';
        errB(3,nx)=T3{nx}(Qx(:,1),Qx(:,2))'*Qw';
    end
    errB=abs(errB-E(:,1));
   
    test(inc2)=min(reshape(err(:,1)./errB,[],1));
    if test(inc2)>best
        fprintf('%i %f %f \n',inc2,test(inc2),conder(inc2)) 
        errs(:,:,inc)=abs(err(:,1)./errB);
        QxM(:,:,inc)=Qx;
        QwM(inc,:)=Qw;
        inc=inc+1;
        best=test(inc2);
    end
    inc2=inc2+1;
end

%Qx Tx=0.0469 Ty=0 y-Ty err=3.5x
% 0.0037    0.7092
% 0.0993    0.5970
% 0.3754    0.9824
% 0.3021    0.4312
% 0.4130    0.0040
% 0.5015    0.9217
% 0.7077    0.7367
% 0.4541    0.9870
% 0.9897    0.1806
% 0.6880    0.6429
% 0.1566    0.4132
% 0.5713    0.6035
% 0.1697    0.1197
% 0.0701    0.2816
% 0.9490    0.6785
% 0.4154    0.9870
% 0.0256    0.3011
% 0.2856    0.2172
% 0.5303    0.0897
% 0.8057    0.7298
% 0.5750    0.5065
%Qw
%[0.0970135821138229,-0.264160266848263,2.56500670267915,-0.331046014174798,-0.0980054414973220,0.308401066616790,-1.40854595203640,2.63501798518489,0.0619952074075802,-0.00654645887268971,0.769261541518412,0.795877905209022,0.264114141373189,-0.398583014379401,-0.0887003651354474,-5.13901753092481,0.143898430481680,-0.199677360854949,0.277857637224571,1.04602233496412,-0.0301841300477840;]