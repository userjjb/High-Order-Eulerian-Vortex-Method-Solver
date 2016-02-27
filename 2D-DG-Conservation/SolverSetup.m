Estreamy= Estreamy; %Make Matlab explicitly aware this is a variable

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
switch KernelType
case 'BS'
    Rsq= sum(bsxfun(@minus,rv_xB,[srcx(:),srcy(:)]).^2,2);
    gkernel_xB= (1/(2*pi))*permute(bsxfun(@minus,srcy(:),rv_xB(1,2,:,:))./Rsq,[1 4 3 2]);
    Rsq= sum(bsxfun(@minus,rv_yB,[srcx(:),srcy(:)]).^2,2);
    gkernel_yB= (1/(2*pi))*permute(bsxfun(@minus,rv_yB(1,1,:,:),srcx(:))./Rsq,[1 4 3 2]);
case 'WL'
    Rsq= sum(bsxfun(@minus,rv_xB,[srcx(:),srcy(:)]).^2,2);
    gkernel_xB= (1/(2*pi))*permute(( bsxfun(@minus,srcy(:),rv_xB(1,2,:,:)).*(Rsq+(2)*del^2) )./(Rsq+del^2).^(2),[1 4 3 2]);
    Rsq= sum(bsxfun(@minus,rv_yB,[srcx(:),srcy(:)]).^2,2);
    gkernel_yB= (1/(2*pi))*permute(( bsxfun(@minus,rv_yB(1,1,:,:),srcx(:)).*(Rsq+(2)*del^2) )./(Rsq+del^2).^(2),[1 4 3 2]);
case 'SG'
    Rsq= sum(bsxfun(@minus,rv_xB,[srcx(:),srcy(:)]).^2,2);
    gkernel_xB= (1/(2*pi))*permute( (bsxfun(@minus,srcy(:),rv_xB(1,2,:,:))./(Rsq)).*(1-(1-Rsq/del^2).*exp(-Rsq/del^2)),[1 4 3 2]);
    Rsq= sum(bsxfun(@minus,rv_yB,[srcx(:),srcy(:)]).^2,2);
    gkernel_yB= (1/(2*pi))*permute( (bsxfun(@minus,rv_yB(1,1,:,:),srcx(:))./(Rsq)).*(1-(1-Rsq/del^2).*exp(-Rsq/del^2)),[1 4 3 2]);
case 'SG6'
    Rsq= sum(bsxfun(@minus,rv_xB,[srcx(:),srcy(:)]).^2,2);
    gkernel_xB= (1/(2*pi))*permute( (bsxfun(@minus,srcy(:),rv_xB(1,2,:,:))./(Rsq)).*(1-(1-2*Rsq/del^2+0.5*Rsq.^2/del^4).*exp(-Rsq/del^2)),[1 4 3 2]);
    Rsq= sum(bsxfun(@minus,rv_yB,[srcx(:),srcy(:)]).^2,2);
    gkernel_yB= (1/(2*pi))*permute( (bsxfun(@minus,rv_yB(1,1,:,:),srcx(:))./(Rsq)).*(1-(1-2*Rsq/del^2+0.5*Rsq.^2/del^4).*exp(-Rsq/del^2)),[1 4 3 2]);
case 'PS'
    Rsq= sum(bsxfun(@minus,rv_xB,[srcx(:),srcy(:)]).^2,2);
    gkernel_xB= (1/(2*pi))*permute( (bsxfun(@minus,srcy(:),rv_xB(1,2,:,:))./(Rsq)).*(1-besselj(0,sqrt(Rsq)/del)),[1 4 3 2]);
    Rsq= sum(bsxfun(@minus,rv_yB,[srcx(:),srcy(:)]).^2,2);
    gkernel_yB= (1/(2*pi))*permute( (bsxfun(@minus,rv_yB(1,1,:,:),srcx(:))./(Rsq)).*(1-besselj(0,sqrt(Rsq)/del)),[1 4 3 2]);
case 'PS2'
    Rsq= sum(bsxfun(@minus,rv_xB,[srcx(:),srcy(:)]).^2,2);
    gkernel_xB= (1/(2*pi))*permute( (bsxfun(@minus,srcy(:),rv_xB(1,2,:,:))./(Rsq)).*(1- (8./(45*Rsq/del^2)).*(4*besselj(2,4*sqrt(Rsq)/del)-5*besselj(2,2*sqrt(Rsq)/del)+besselj(2,sqrt(Rsq)/del))),[1 4 3 2]);
    Rsq= sum(bsxfun(@minus,rv_yB,[srcx(:),srcy(:)]).^2,2);
    gkernel_yB= (1/(2*pi))*permute( (bsxfun(@minus,rv_yB(1,1,:,:),srcx(:))./(Rsq)).*(1- (8./(45*Rsq/del^2)).*(4*besselj(2,4*sqrt(Rsq)/del)-5*besselj(2,2*sqrt(Rsq)/del)+besselj(2,sqrt(Rsq)/del))),[1 4 3 2]);
otherwise
    gkernel_xB= (1/(2*pi))*permute(bsxfun(@minus,srcy(:),rv_xB(1,2,:,:))./(sum(bsxfun(@minus,rv_xB,[srcx(:),srcy(:)]).^2,2)+del^2),[1 4 3 2]);
    gkernel_yB= (1/(2*pi))*permute(bsxfun(@minus,rv_yB(1,1,:,:),srcx(:))./(sum(bsxfun(@minus,rv_yB,[srcx(:),srcy(:)]).^2,2)+del^2),[1 4 3 2]);
end

%% Near kernel and local params
%Function for selecting nearby elements in range, local stream coords
LEnum= reshape(1:(2*NearRange+1)^2,2*NearRange+1,[]);
Local= @(source) LEnum( (max(Enumy(source)-NearRange,1):min(Enumy(source)+NearRange,K(2)))-Enumy(source)+NearRange+1 ,...
    (max(Enumx(source)-NearRange,1):min(Enumx(source)+NearRange,K(1)))-Enumx(source)+NearRange+1 );
Lstreamx= reshape(reshape(1:Np*(2*NearRange+1)^2,2*NearRange+1,[])',Np,[]);
Lstreamy= reshape(permute(reshape(reshape(1:Np*(2*NearRange+1)^2,2*NearRange+1,[])',Np,2*NearRange+1,[]),[1 3 2]),Np,[]);

%Function for selecting nearby elements in range, global stream coords
Near=@(source) Enum(max(Enumy(source)-NearRange,1):min(Enumy(source)+NearRange,end),...
    max(Enumx(source)-NearRange,1):min(Enumx(source)+NearRange,end));

Nsx=zeros(Np*(2*NearRange+1)^2,Enum(end,end)); Nsy=Nsx; Lsx=Nsx; Lsy=Lsx;
for Src=1:Enum(end,end)
    numS(Src)=Np*numel(Near(Src));
    Nsx(1:numS(Src),Src)= sort(reshape(Estreamx(:, Near(Src)),[],1));
    Nsy(1:numS(Src),Src)= sort(reshape(Estreamy(:, Near(Src)),[],1));
    Lsx(1:numS(Src),Src)= sort(reshape(Lstreamx(:, Local(Src)),[],1));
    Lsy(1:numS(Src),Src)= sort(reshape(Lstreamy(:, Local(Src)),[],1));
end
%Note: Nrv y-coords for x-streams are the same as x-coords for y-streams and vice versa
[Nrv_y,Nrv_x]= meshgrid(reshape(bsxfun(@plus,Qx/2,-NearRange:NearRange)*delX,1,[]),...
    reshape(bsxfun(@plus,Qx2(2:end-1)/2,-NearRange:NearRange)*delX,1,[]));
[srcx,srcy]= meshgrid(Qx*(delX/2),Qx*(delX/2));

Rsqx= bsxfun(@minus,Nrv_x(:)',srcx(:)).^2 + bsxfun(@minus,Nrv_y(:)',srcy(:)).^2;
Rsqy= bsxfun(@minus,Nrv_y(:)',srcx(:)).^2 + bsxfun(@minus,Nrv_x(:)',srcy(:)).^2;
Mx= bsxfun(@minus,srcy(:),Nrv_y(:)')/(2*pi);
My= bsxfun(@minus,Nrv_y(:)',srcx(:))/(2*pi);
switch KernelType
case 'BS'
    kernel_x= Mx./Rsqx;
    kernel_y= My./Rsqy;
case 'WL'
    kernel_x= ( Mx.*(Rsqx+2*del^2) )./(Rsqx+del^2).^(2);
    kernel_y= ( My.*(Rsqy+2*del^2) )./(Rsqy+del^2).^(2);
case 'SG'
    kernel_x= (Mx./Rsqx).*(1-(1-Rsqx/del^2).*exp(-Rsqx/del^2));
    kernel_y= (My./Rsqy).*(1-(1-Rsqy/del^2).*exp(-Rsqy/del^2));
case 'SG6'
    kernel_x= (Mx./Rsqx).*(1-(1-2*Rsqx/del^2+0.5*Rsqx.^2/del^4).*exp(-Rsqx/del^2));
    kernel_y= (My./Rsqy).*(1-(1-2*Rsqy/del^2+0.5*Rsqy.^2/del^4).*exp(-Rsqy/del^2));
case 'PS'
    kernel_x= (Mx./Rsqx).*(1-besselj(0,sqrt(Rsqx)/del));
    kernel_y= (My./Rsqy).*(1-besselj(0,sqrt(Rsqy)/del));
case 'PS2'
    kernel_x= (Mx./Rsqx).*(1- (8./(45*Rsqx/del^2)).*(4*besselj(2,4*sqrt(Rsqx)/del)-5*besselj(2,2*sqrt(Rsqx)/del)+besselj(2,sqrt(Rsqx)/del)));
    kernel_y= (My./Rsqy).*(1- (8./(45*Rsqy/del^2)).*(4*besselj(2,4*sqrt(Rsqy)/del)-5*besselj(2,2*sqrt(Rsqy)/del)+besselj(2,sqrt(Rsqy)/del)));
otherwise %RM
    kernel_x= Mx./(Rsqx+del^2);
    kernel_y= My./(Rsqy+del^2);
end
%clearvars Nrv_x Nrv_y srcx srcy Local LEnum Lstreamx Lstreamy Rsq Rsqx Rsqy
%% Source and neighbor elements modified kernel
%The singular (for self), or nearly-singular (for neighboring elements) BS kernel leads to
%unallowable quadrature error. We can generate modified kernel values that when used lead to correct
%integration of these elements. They are modified here.
if strcmp(KernelType,'BS')
load(ModKernelFile); %Load Iwm
%This file is expected to contain 'Iwm' a (Mp,Np,Mp,Np,3,3,2) size matrix. The value ordering is 
%(J,I,Tj,Ti,Ex,Ey,X/Y) where J,I is the associated source vorticity interp node; Tj,Ti is the target
%velocity stream node in the target element; Ex,Ey is the target element number with 2,2 being the
%self element and 1,1 is the element with both x and y coordinates less than the self element; the
%7th dimension has x-velocity kernels in the first entry and y-velocity kernel values in the second

NBnum= bsxfun(@plus, [1:3*(Mp-2)]', ((NearRange-1)*(Mp-2)*((2*NearRange+1)*Np+1)) + ((2*NearRange+1)*(Mp-2)*[0:3*Np-1]));
kernMap= reshape(Qw'*Qw,[],1)*((B(2)-B(1))/(K(1)*4));
kernel_x(:,NBnum(:))= bsxfun(@rdivide,reshape(permute(Iwm(:,:,:,2:end-1,:,:,1),[1 2 4 5 3 6 7]),[Np^2 3*Np*3*(Mp-2)]),kernMap);
kernel_y(:,NBnum(:))= bsxfun(@rdivide,reshape(permute(Iwm(:,:,2:end-1,:,:,:,2),[1 2 3 6 4 5 7]),[Np^2 3*Np*3*(Mp-2)]),kernMap);
end
%Near kernel modification

%%

%Outer product of vorticity quadrature weights for pre-multiplication,
%including Jacobian
QwPre=(delX/2)^2*reshape(Qw'*Qw,[],1);
k2= zeros(size(wx)); %LSERK stage state
w_tot= abs(permute(mtimesx(reshape(permute(reshape(wy,Np,K(2),Np,K(1)),[1 3 2 4]),1,Np^2,K(2)*K(1)),QwPre),[3 1 2]));
mask= find(w_tot>w_thresh);
setup= [sum(w_tot),N,M,del,delt,EndTime,K(1),K(2),B,TestCases,NearRange,0];
zmax= 1.5*max(max(w)); zmin=1.5*min(min(w));
itt= 1; BackupSave= 1800;
StepNum= uint64(0);