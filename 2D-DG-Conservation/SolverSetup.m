Estreamy= Estreamy; %Make Matlab explicitly aware this is a variable

%Function for selecting nearby elements in range
Near=@(source) Enum(max(Enumy(source)-NearRange,1):min(Enumy(source)+NearRange,end),...
    max(Enumx(source)-NearRange,1):min(Enumx(source)+NearRange,end));
%Setup:
%Weighting basis with lumped outside terms to save on calculation for
%numerical total nodal surface flux term: \hat{f}_R L_j(x_R)-\hat{f}_L L_j(x_L)
%Distributed: inv Mass matrix (1/QwG), inv Jacobian (2/delX), and 1/2 coeff
%from num flux eval
LlM=  permute( Ll.*(1./(delX*Qw')) , [2 3 4 1]); %Left elem boundary eval
LrM=  permute( Lr.*(1./(delX*Qw')) , [2 3 4 1]); %Right elem boundary eval

wx=   reshape(w',Np,1,[]);      %Reshape vorticity for mtimesx_x bsx
wy=   reshape(w,Np,1,[]);       %Reshape vorticity for mtimesx_y bsx

%Intialize scalar kernel values, a specific element stream's kernel should be a
%matrix of pointwise values: each column corresponds to the matching point
%in the stream, each element in the column matches the scalar value of the
%kernel for that point as one moves col-wise within the source
%element to match w_elem's ordering.
srcx=wxm(Nnumy(:,1,Estreamy(:,end) )); %Generalized source location
srcy=wym(Nnumy(:,1,Estreamy(:,end) ));
%Calculate generalized kernel for boundary velocity points
gkernel_xB= (1/(4*pi))*permute(bsxfun(@minus,srcy(:),rv_xB(1,2,:,:))./(sum(bsxfun(@minus,rv_xB,[srcx(:),srcy(:)]).^2,2)+del^2).^(3/2),[1 4 3 2]);
gkernel_yB= (1/(4*pi))*permute(bsxfun(@minus,rv_yB(1,1,:,:),srcx(:))./(sum(bsxfun(@minus,rv_yB,[srcx(:),srcy(:)]).^2,2)+del^2).^(3/2),[1 4 3 2]);

rv_xS=zeros(1,2,Mp-2,Np*(2*NearRange+1)^2,Enum(end,end)); rv_yS=rv_xS;
Nsx=zeros(Np*(2*NearRange+1)^2,Enum(end,end)); Nsy=Nsx;
for Src=1:Enum(end,end)
    numS(Src)=Np*numel(Near(Src));
    Nsx(1:numS(Src),Src)= sort(reshape(Estreamx(:, Near(Src)),[],1));
    Nsy(1:numS(Src),Src)= sort(reshape(Estreamy(:, Near(Src)),[],1));
    rv_xS(:,:,:,1:numS(Src),Src)= rv_x(:,:,:,Nsx(1:numS(Src),Src) );
    rv_yS(:,:,:,1:numS(Src),Src)= rv_y(:,:,:,Nsy(1:numS(Src),Src) );
end

srcx=wxm(reshape(Nnumy(:,1,Estreamy),Np^2,1,1,1,[]));
srcy=wym(reshape(Nnumy(:,1,Estreamy),Np^2,1,1,1,[]));
kernel_x= (1/(4*pi))*permute(bsxfun(@minus,srcy,rv_xS(1,2,:,:,:))./(sum(bsxfun(@minus,rv_xS,[srcx,srcy]).^2,2)+del^2).^(3/2),[1 3 4 5 2]);
kernel_y= (1/(4*pi))*permute(bsxfun(@minus,rv_yS(1,1,:,:,:),srcx)./(sum(bsxfun(@minus,rv_yS,[srcx,srcy]).^2,2)+del^2).^(3/2),[1 3 4 5 2]);
clearvars srcx srcy %rv_xS rv_yS

%Outer product of vorticity quadrature weights for pre-multiplication,
%including Jacobian
QwPre=(delX/2)^2*reshape(Qw'*Qw,1,[]);
k2= zeros(size(wx)); %LSERK stage state
w_tot=abs(permute(mtimesx(reshape(permute(reshape(wy,Np,K(2),Np,K(1)),[1 3 2 4]),1,Np^2,K(2)*K(1)),QwPre'),[3 1 2]));
mask=find(w_tot>w_thresh);
setup=[sum(w_tot),N,M,del,delt,EndTime,K(1),K(2),B,TestCases];
zmax=1.5*max(max(w)); zmin=1.5*min(min(w));
itt=0;