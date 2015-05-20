clear all
clc

load('6_12GDdt2.mat')
w_tot0=setup(1); N=setup(2); M=setup(3); del=setup(4); delt2=setup(5); EndTime=setup(6); K=setup(7:8); B=setup(9:12);
run('CalcedParams')
Np2=Np; QxA=Qx; Estreamx2=Estreamx; Enumx2=Enumx; Enumy2=Enumy; Enum2=Enum;
wxA=permute(reshape(permute(wxt(:,1,Estreamx,:),[1 3 2 4]),Np,Np,K(2),K(1),[]),[2 1 3 4 5]);

load('9_36GD.mat')
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

Id=12;
Ip= 1/(2*Id) : 1/Id : 1-1/(2*Id);

Id2=36;
Ip2= 1/(2*Id2) : 1/Id2 : 1-1/(2*Id2);

interp_g=@(x,y,G) Lag(x,1:Np)'*(G*Lag(y,1:Np));
Lx= Lag(Ip,1:Np)';
Ly= Lag(Ip,1:Np);
Lx2= Lag2(Ip2,1:Np2)';
Ly2= Lag2(Ip2,1:Np2);

freq= lcm(delt2*100,delt*100)/(delt*100);
freq2= lcm(delt2*100,delt*100)/(delt2*100);
tratio= freq2/freq;
endt=size(wxE,5);
endt2=size(wxA,5);

h=waitbar(0);
it=1;
tic
for t=1:freq:endt
    wxIA= permute(mtimesx(Lx2,mtimesx( wxA(:,:,:,:,1+(t-1)*tratio),Ly2 )),[1 3 2 4]);
    wxIE= permute(mtimesx(Lx, mtimesx( wxE(:,:,:,:,t),Ly )),[1 3 2 4]);

    L2(it)= sqrt((delX/Id)^2*sum((wxIE(:)-wxIA(:)).^2));
    tt(it)= delt*(t-1);
    waitbar(t/endt,h,sprintf('Estimated time: %i',floor((endt-t)*toc/t)))
    it= it+1;
end
close(h)
plot(tt,L2,'k')
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