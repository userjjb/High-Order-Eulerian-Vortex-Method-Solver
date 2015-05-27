clear all
clc

load('6_12_100_Stage.mat')
w_tot0=setup(1); N=setup(2); M=setup(3); del=setup(4); delt=setup(5); EndTime=setup(6); K=setup(7:8); B=setup(9:12);
run('CalcedParams')
Estreamx2=Estreamx; Enumx2=Enumx; Enumy2=Enumy; Enum2=Enum;
wxA=permute(reshape(permute(wxt(:,1,Estreamx,:),[1 3 2 4]),Np,Np,prod(K),[]),[2 1 3 4]);

load('6_48_100_Stage.mat')
w_tot0=setup(1); N=setup(2); M=setup(3); del=setup(4); delt=setup(5); EndTime=setup(6); K=setup(7:8); B=setup(9:12);
run('CalcedParams')
wxE=permute(reshape(permute(wxt(:,1,Estreamx,:),[1 3 2 4]),Np,Np,prod(K),[]),[2 1 3 4]);

clearvars wxt

nd=(Qx+1)/2;
nn=elim(nd',nd',[1 3 2]);
Lag= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,nn(nv,:,:)),bsxfun(@minus,nd(nv),nn(nv,:,:))),3);

Id=16;
Ip= 1/(2*Id) : 1/Id : 1-1/(2*Id);

Id2=64;
Ip2= 1/(2*Id2) : 1/Id2 : 1-1/(2*Id2);

interp_g=@(x,y,G) Lag(x,1:Np)'*(G*Lag(y,1:Np));
Lx= Lag(Ip,1:Np)';
Ly= Lag(Ip,1:Np);
Lx2= Lag(Ip2,1:Np)';
Ly2= Lag(Ip2,1:Np);

h=waitbar(0);
tic
for t=1:size(wxE,4)
    wxIE= mtimesx(Lx, mtimesx( wxE(:,:,:,t),Ly ) );
    wxIA= mtimesx(Lx2, mtimesx( wxA(:,:,:,t),Ly2 ) );
    
    for i=1:48
    for j=1:48
        wxIEM(Id*(i-1)+(1:Id),Id*(j-1)+(1:Id))=wxIE(:,:,i+(j-1)*48);
    end
    end
    for i=1:12
    for j=1:12
        wxIAM(Id2*(i-1)+(1:Id2),Id2*(j-1)+(1:Id2))=wxIA(:,:,i+(j-1)*12);
    end
    end
    L2(t)= sqrt((delX/Id)^2*sum((wxIEM(:)-wxIAM(:)).^2));
    waitbar(t/size(wxE,4),h,sprintf('Estimated time: %i',floor((size(wxE,4)-t)*toc/t)))
end
close(h)
plot(L2,'g')
hold on


% IpD=bsxfun(@plus,Ip'*delX,Ex(1:end-1));
% [xx,yy]=meshgrid(IpD(:),IpD(:));
% 
% contour(xx,yy,wxI,[-.946:.07:.45],'linewidth',.3);
% axis equal; axis(B);