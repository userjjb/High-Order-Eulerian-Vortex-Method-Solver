clear all
clc

load('6_36G12_del.mat')
w_tot0=setup(1); N=setup(2); M=setup(3); del=setup(4); delt=setup(5); EndTime=setup(6); K=setup(7:8); B=setup(9:12);
run('CalcedParams')
Np2=Np; QxA=Qx; Estreamx2=Estreamx; Enumx2=Enumx; Enumy2=Enumy; Enum2=Enum;
wxA=permute(reshape(permute(wxt(:,1,Estreamx,:),[1 3 2 4]),Np,Np,K(2),K(1),[]),[2 1 3 4 5]);

load('6_48_100_Stage_NR16.mat')
w_tot0=setup(1); N=setup(2); M=setup(3); del=setup(4); delt=setup(5); EndTime=setup(6); K=setup(7:8); B=setup(9:12);
run('CalcedParams')
wxE=permute(reshape(permute(wxt(:,1,Estreamx,:),[1 3 2 4]),Np,Np,K(2),K(1),[]),[2 1 3 4 5]);

clearvars wxt

nd=(Qx+1)/2;
nn=elim(nd',nd',[1 3 2]);
Lag= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,nn(nv,:,:)),bsxfun(@minus,nd(nv),nn(nv,:,:))),3);
nd2=(QxA+1)/2;
nn2=elim(nd2',nd2',[1 3 2]);
Lag2= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,nn2(nv,:,:)),bsxfun(@minus,nd2(nv),nn2(nv,:,:))),3);

Id=15;
Ip= 1/(2*Id) : 1/Id : 1-1/(2*Id);

Id2=20;
Ip2= 1/(2*Id2) : 1/Id2 : 1-1/(2*Id2);

interp_g=@(x,y,G) Lag(x,1:Np)'*(G*Lag(y,1:Np));
Lx= Lag(Ip,1:Np)';
Ly= Lag(Ip,1:Np);
Lx2= Lag2(Ip2,1:Np2)';
Ly2= Lag2(Ip2,1:Np2);

h=waitbar(0);
endt=size(wxE,5);
wxIE= permute(mtimesx(Lx, mtimesx( wxE(:,:,:,:,285),Ly )),[1 3 2 4]);
tic
for t=1:endt
    wxIA= permute(mtimesx(Lx2,mtimesx( wxA(:,:,:,:,t),Ly2 )),[1 3 2 4]);
    

    L2(t)= sqrt((delX/Id)^2*sum((wxIE(:)-wxIA(:)).^2));
    waitbar(t/endt,h,sprintf('Estimated time: %i',floor((endt-t)*toc/t)))
end
close(h)
plot(L2,'c')
hold on

% endt=1+round((size(wxE,5)-1)/5);
% tic
% for t=1:endt
%     wxIE= permute(mtimesx(Lx,mtimesx( wxE(:,:,:,:,5*(t-1)+1),Ly )),[1 3 2 4]);
%     wxIA= permute(mtimesx(Lx2, mtimesx( wxA(:,:,:,:,t),Ly2 )),[1 3 2 4]);
% 
%     L2(t)= sqrt((delX/Id)^2*sum((wxIE(:)-wxIA(:)).^2));
%     waitbar(t/endt,h,sprintf('Estimated time: %i',floor((endt-t)*toc/t)))
% end

% IpD=bsxfun(@plus,Ip'*delX,Ex(1:end-1));
% [xx,yy]=meshgrid(IpD(:),IpD(:));
% 
% contour(xx,yy,wxI,[-.946:.07:.45],'linewidth',.3);
% axis equal; axis(B);