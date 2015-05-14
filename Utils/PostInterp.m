clear all
clc

load('6_48_100_Stage.mat')
w_tot0=setup(1); N=setup(2); M=setup(3); del=setup(4); delt=setup(5); EndTime=setup(6); K=setup(7:8); B=setup(9:12);
run('CalcedParams')
Estreamx2=Estreamx; Enumx2=Enumx; Enumy2=Enumy; Enum2=Enum;
wxA=permute(reshape(permute(wxt(:,1,Estreamx,:),[1 3 2 4]),Np,Np,prod(K),[]),[2 1 3 4]);

load('6-7_48_100_Stage.mat')
w_tot0=setup(1); N=setup(2); M=setup(3); del=setup(4); delt=setup(5); EndTime=setup(6); K=setup(7:8); B=setup(9:12);
run('CalcedParams')
wxE=permute(reshape(permute(wxt(:,1,Estreamx,:),[1 3 2 4]),Np,Np,prod(K),[]),[2 1 3 4]);

clearvars wxt

nd=(Qx+1)/2;
nn=elim(nd',nd',[1 3 2]);
Lag= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,nn(nv,:,:)),bsxfun(@minus,nd(nv),nn(nv,:,:))),3);

Ip= 0.0625/2:.0625:1-0.0625/2;
NI= numel(Ip);

frac=1;
Ip2= 0.0625/(2*frac):.0625/frac:1-0.0625/(2*frac);
NI2= numel(Ip2);

interp_g=@(x,y,G) Lag(x,1:Np)'*(G*Lag(y,1:Np));
Lx= Lag(Ip,1:Np)';
Ly= Lag(Ip,1:Np);
Lx2= Lag(Ip2,1:Np)';
Ly2= Lag(Ip2,1:Np);

for t=1:size(wxE,4)
    wxIE= mtimesx(Lx, mtimesx( wxE(:,:,:,t),Ly ) );
    wxIA= mtimesx(Lx2, mtimesx( wxA(:,:,:,t),Ly2 ) );
    L2(t)= sqrt((delX/16)^2*sum((wxIE(:)-wxIA(:)).^2));    
end
plot(L2,'g')
hold on

% IpD=bsxfun(@plus,Ip'*delX,Ex(1:end-1));
% [xx,yy]=meshgrid(IpD(:),IpD(:));
% 
% contour(xx,yy,wxI,[-.946:.07:.45],'linewidth',.3);
% axis equal; axis(B);