clear all
clc

%-------------
load('6_12_100_Stage4.mat')
wxA=wxt;
w_tot0=setup(1); N=setup(2); M=setup(3); del=setup(4); delt=setup(5); EndTime=setup(6); K=setup(7:8); B=setup(9:12);
run('CalcedParams')
Estreamx2=Estreamx; Enumx2=Enumx; Enumy2=Enumy; Enum2=Enum;

Ip2= 0.0625/8:.0625/4:1-0.0625/8;
NI2= numel(Ip2);

Wx2 = LagBaryWeight(Qx);
[xx,yy]=meshgrid(Ip2,Ip2);
Reval=[xx(:),yy(:)];
nR= size(Reval,1); %Number of eval points
nd= (Qx+1)/2;

%We essentially would like to calculate lX* sum: Wi*zi/deni as a 1D element-wise example
denx= bsxfun(@minus,Reval(:,1),nd'); %Basis specific denominator that "fixes" lX by dividing it to exclude the unwanted (x-x_i) for basis 'ij'
deny= bsxfun(@minus,Reval(:,2),nd');
lX2= prod(denx,2); %The numerator in the sum
lY2= prod(deny,2);
%den=|-------------------denx-------|.*|----deny------|
den2= reshape(repmat(denx,Np,1),nR,Np^2).*repmat(deny,1,Np); %Forms total den, we saved denx and deny only so we could calc lX and lY

Wy2=repmat(Wx2,1,Np);
Wx2=reshape(repmat(Wx2,Np,1),1,Np^2); %Form the vectors to calc the tensor product of Wi*Wj

%-------------
load('6_48_100_Stage_NR16.mat')
wxE=wxt;
w_tot0=setup(1); N=setup(2); M=setup(3); del=setup(4); delt=setup(5); EndTime=setup(6); K=setup(7:8); B=setup(9:12);
run('CalcedParams')
clearvars wxt

Ip= 0.0625/2:.0625:1-0.0625/2;
NI= numel(Ip);
Ip3= 0.0625/4:.0625/2:1-0.0625/4;
NI3= numel(Ip3);

Wx = LagBaryWeight(Qx);
[xx,yy]=meshgrid(Ip,Ip);
Reval=[xx(:),yy(:)];
nR= size(Reval,1); %Number of eval points
nd= (Qx+1)/2;

%We essentially would like to calculate lX* sum: Wi*zi/deni as a 1D element-wise example
denx= bsxfun(@minus,Reval(:,1),nd'); %Basis specific denominator that "fixes" lX by dividing it to exclude the unwanted (x-x_i) for basis 'ij'
deny= bsxfun(@minus,Reval(:,2),nd');
lX= prod(denx,2); %The numerator in the sum
lY= prod(deny,2);
%den=|-------------------denx-------|.*|----deny------|
den= reshape(repmat(denx,Np,1),nR,Np^2).*repmat(deny,1,Np); %Forms total den, we saved denx and deny only so we could calc lX and lY

Wy=repmat(Wx,1,Np);
Wx=reshape(repmat(Wx,Np,1),1,Np^2); %Form the vectors to calc the tensor product of Wi*Wj

%------------------


h=waitbar(0);
tic
for t=1:size(wxE,4)
    for Src=1:Enum(end,end)
        wxIE( NI*(Enumy(Src)-1)+(1:NI), NI*(Enumx(Src)-1)+(1:NI) )= reshape( lX.*(bsxfun(@times,Wx,bsxfun(@rdivide,Wy,den))*reshape(permute(wxE(:,1,Estreamx(:,Src),t),[3 1 2]),[],1)).*lY ,NI,NI);
    end
    for Src=1:Enum2(end,end)
        wxIA( NI2*(Enumy2(Src)-1)+(1:NI2), NI2*(Enumx2(Src)-1)+(1:NI2) )= reshape( lX2.*(bsxfun(@times,Wx2,bsxfun(@rdivide,Wy2,den2))*reshape(permute(wxA(:,1,Estreamx2(:,Src),t),[3 1 2]),[],1)).*lY2 ,NI2,NI2);
    end

    L2(t)= sqrt((delX/16)^2*sum(sum((wxIE-wxIA).^2)));    
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