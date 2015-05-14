clear all
clc

load('6-7_48_100_Stage.mat')
wxA=wxt;
w_tot0=setup(1); N=setup(2); M=setup(3); del=setup(4); delt=setup(5); EndTime=setup(6); K=setup(7:8); B=setup(9:12);
run('CalcedParams')
Estreamx2=Estreamx; Enumx2=Enumx; Enumy2=Enumy; Enum2=Enum;

load('6_48_100_Stage.mat')
wxE=wxt;
w_tot0=setup(1); N=setup(2); M=setup(3); del=setup(4); delt=setup(5); EndTime=setup(6); K=setup(7:8); B=setup(9:12);
run('CalcedParams')

clearvars wxt

nd=(Qx+1)/2;
nn=elim(nd',nd',[1 3 2]);
Lag= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,nn(nv,:,:)),bsxfun(@minus,nd(nv),nn(nv,:,:))),3);

Ip= 0.0625/2:.0625:1-0.0625/2;
NI= numel(Ip);
Ip2= 0.0625/2:.0625/1:1-0.0625/2;
NI2= numel(Ip2);
Ip3= 0.0625/4:.0625/2:1-0.0625/4;
NI3= numel(Ip3);

interp_g=@(x,y,G) Lag(x,1:Np)'*(G*Lag(y,1:Np));

h=waitbar(0);
tic
for t=1:size(wxE,4)
    for Src=1:Enum(end,end)
        wxIE( NI*(Enumy(Src)-1)+(1:NI), NI*(Enumx(Src)-1)+(1:NI) )= interp_g(Ip,Ip, permute(wxE(:,1,Estreamx(:,Src),t),[3 1 2]) );
    end
    for Src=1:Enum2(end,end)
        wxIA( NI2*(Enumy2(Src)-1)+(1:NI2), NI2*(Enumx2(Src)-1)+(1:NI2) )= interp_g(Ip2,Ip2, permute(wxA(:,1,Estreamx2(:,Src),t),[3 1 2]) );
    end

    L2(t)= sqrt((delX/16)^2*sum(sum((wxIE-wxIA).^2)));    
    waitbar(t/size(wxE,4),h,sprintf('Estimated time: %i',floor((size(wxE,4)-t)*toc/t)))
end
close(h)
plot(L2,'r')
hold on

% IpD=bsxfun(@plus,Ip'*delX,Ex(1:end-1));
% [xx,yy]=meshgrid(IpD(:),IpD(:));
% 
% contour(xx,yy,wxI,[-.946:.07:.45],'linewidth',.3);
% axis equal; axis(B);