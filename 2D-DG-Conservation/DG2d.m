%Some terminology: the structured tensor grid means we can completely
%separate x and y components and discretize independently. They only couple
%when we calculate dw_p/dt at a point 'p' from the sum of each directions
%contribution, so dw_p/dt = dw_px/dt + dw_py/dt
%The set of colinear points along a coordinate direction will be called a
%"stream"

close all
clear all
clc
%Solver parameters
alpha= 1;                           %Numerical flux param (1 upwind,0 CD)
N= 5;                               %Local vorticity poly order
M= 5;                               %Local velocity poly order
[RKa,RKb,RKc,nS]= LSRKcoeffs('NRK14C');
w_thresh=1E-7;
del=2*0.15^2;
delt= 0.025;
skip= 1;
endtime=28;
DGmask='full';
BCtype= 'NoInflow';
TestCase=2;
%---Global domain initialization (parameters)------------------------------
B= 3.5*[-1.25 1 -1.25 1];                     %left, right, bottom, top
K= [32 32];                         %Num elements along x,y
Ex= linspace(B(1),B(2),K(1)+1);     %Elem edges left-right
Ey= linspace(B(3),B(4),K(2)+1);     %Elem edges bottom-top

%---Global domain initialization (calculated)------------------------------
Np= N+1;                            %Number of interpolation points to achieve N
Mp= M+1;                            %Number of interpolation points to achieve M
delX= (B(2)-B(1))/K(1);             %Element size
    %For simplicity we require square elements
    assert(abs(delX-(B(4)-B(3))/K(2))<eps(delX),'Elements must be square')
    
[Qx, Qw]= GLquad(Np);               %Gauss-Legendre quadrature nodes/weights
[Qx2, Qw2]= LGLquad(Mp);            %Lobatto-Gauss-Legendre
[QwS, Ll, Lr]= StiffnessQuadModWeights(Qx,Qx2); %Modified quadrature weights for the stiffness matrix calculation, and vorticity interpolation evals
QwSM=bsxfun(@times,QwS,2./(delX*permute(Qw,[1 4 3 2]))); %Distributed Mass matrix and Jacobian
[QwSlow,~,~]= StiffnessQuadModWeights(Qx,[-1;1]);
QwSMlow=bsxfun(@times,QwSlow,2./(delX*permute(Qw,[1 4 3 2])));
    
w= zeros(Np*K(2),Np*K(1));          %Global vorticity at element interp points
%NOTE: Be careful to notice the use of Mp vs Np here
%For velocity grids that are non-square tensor products (ie. Mp=!Np) v_x is found
%along x-streams and v_y is found along y-streams, but v_x eval points are
%not generally coincident with v_y eval points. If one needed *both* v_x
%and v_y for a given stream this would be a problem, one would need double
%the velocity evaluations (w/o the benefit of coincident eval points).
%However for stiffness calcs one only needs the velocity component parallel
%to the stream direction. The Biot-Savart integral calculates velocity
%componentwise, so we can take advantage of this and avoid half of the
%velocity evals.
v_xI=zeros(1,Mp-2,K(1)*Np*K(2));
v_yI=zeros(1,Mp-2,K(2)*Np*K(1));
v_xB=zeros(Np*K(2),K(1)+1);
v_yB=zeros(K(2)+1,Np*K(1));

%---Node numbering---------------------------------------------------------
%Element numbering
[Enumx,Enumy]=meshgrid(1:K(1),1:K(2));
%Stream global numbering as well as stream to left or right for periodic
%BCs For instance,pass x_km1 as an index to get the stream left of current
Eord_x= reshape(1:K(1)*K(2)*Np,K(1),K(2)*Np)';
x_km1= reshape(circshift(Eord_x,[0 1])',[],1);
x_kp1= reshape(circshift(Eord_x,[0 -1])',[],1);
Eord_y= reshape(1:K(1)*K(2)*Np,K(2),K(1)*Np);
y_km1= reshape(circshift(Eord_y,[1 0]),[],1);
y_kp1= reshape(circshift(Eord_y,[-1 0]),[],1);
%Boundary velocity node numbering
EBx=reshape(1:Np*K(2)*(K(1)+1),Np*K(2),(K(1)+1));
EBl= reshape(EBx(:,1:end-1)',1,1,[]);
EBr= reshape(EBx(:,2:end)',1,1,[]);
EBy=reshape(1:Np*K(1)*(K(2)+1),(K(2)+1),Np*K(1));
EBb= reshape(EBy(1:end-1,:),1,1,[]);
EBt= reshape(EBy(2:end,:),1,1,[]);
%Stream/element associativity (stream in:,element) col-wise
Estreamx= reshape(Eord_x,Np,[]);
Estreamy= reshape(permute(reshape(Eord_y',Np,K(1),[]),[1 3 2]),Np,[]);
%Vorticity node/stream associativity
Nnum=reshape(1:numel(w),size(w));
Nnumx=reshape(Nnum',Np,1,[]);
Nnumy=reshape(Nnum,Np,1,[]);

%---Node coordinates-------------------------------------------------------
%Interp/quad node locations
x_w= reshape(bsxfun(@plus,repmat((Qx+1)*(delX/2),1,K(1)),Ex(1:end-1)),[],1);
y_w= reshape(bsxfun(@plus,repmat((Qx+1)*(delX/2),1,K(2)),Ey(1:end-1)),[],1);
%Elementwise velocity interp/quad node locations, reflected about axes
[t2,t1]= meshgrid([y_w;(B(4)-B(3))+y_w],[Ex, (B(2)-B(1))+Ex(2:end)]);
rv_xB= reshape([t1(:),t2(:)]',1,2,2*K(1)+1,[]);
[t1,t2]= meshgrid([x_w;(B(2)-B(1))+x_w],[Ey,(B(4)-B(3))+Ey(2:end)]);
rv_yB= reshape([t1(:),t2(:)]',1,2,2*K(2)+1,[]);

x_v=reshape(bsxfun(@plus,(Qx2(2:end-1)+1)*(delX/2),Ex(1:end-1)),1,[]);
y_v=reshape(bsxfun(@plus,(Qx2(2:end-1)+1)*(delX/2),Ey(1:end-1)),1,[]);
[t2,t1]= meshgrid([y_w;(B(4)-B(3))+y_w],[x_v, (B(2)-B(1))+x_v]);
rv_x= reshape([t1(:),t2(:)]',1,2,2*K(1)*(Mp-2),[]);
[t1,t2]= meshgrid([x_w;(B(2)-B(1))+x_w],[y_v, (B(2)-B(1))+y_v]);
rv_y= reshape([t1(:),t2(:)]',1,2,2*K(2)*(Mp-2),[]);
%Note that the actual velocity node locations for x_streams are (x_v,y_w) and
%for y_streams are (y_v,x_w). The velocity grid is NOT a square tensor product, but
%rather velocity "streams" are coincident with vorticity streams.
%Additionally, the boundary velocity locations occur at (Ex,y_w) and
%(x_w,Ey) analogous to the interpolation points
[wxm, wym]= meshgrid(x_w,y_w);

%For use later with discrete norm calcs
h_x=(diff([B(1)+(B(1)-x_w(1));x_w]) + diff([x_w; B(2)+(B(2)-x_w(end))]))/2;
h_y=(diff([B(3)+(B(3)-y_w(1));y_w]) + diff([y_w; B(4)+(B(4)-y_w(end))]))/2;
[h_x, h_y]= meshgrid(h_x,h_y);
norm_h=h_x.*h_y;

%2D mollifier with center dx,dy and range over the ellipse with axes a,b
%Not currently used, but may prove useful for IC funs that have 
%discontinuous derivatives at truncation edges for smoothing to avoid the
%interpolation from going haywire
moll=@(x,y,dx,dy,a,b) heaviside(1-(((x-dx)/a).^2+((y-dy)/b).^2)).*exp(1+(-1./(1-(((x-dx)/a).^2+((y-dy)/b).^2).^2)));
%Intial conditions specification-------------------------------------------
ICfuns={}; %Create a cell list of functions that define the ICs
    %Gaussian 1
Ga1=    0.012;
Gb1=    0.025;
Gdx1=   0.12; 
Gdy1=   0.12;
GA1=    .5;
ICfuns{end+1}=@(x,y) GA1*exp(-(((x)-Gdx1).^2/Ga1+(y-Gdy1).^2/Gb1));%Center is at (dx,dy)
    %Gaussian 2
Ga2=    0.012;
Gb2=    0.025;
Gdx2=   -0.12; 
Gdy2=   -0.12;
GA2=    .5;
ICfuns{end+1}=@(x,y) GA2*exp(-(((x)-Gdx2).^2/Ga2+(y-Gdy2).^2/Gb2));
    %Fun 3
Pa=1;
Pb=0.5;
PR=.7;
ICfuns{end+1}=@(x,y) (1/PR^14)*(PR^2-min((x/Pa).^2+(y/Pb).^2,PR^2)).^7;
    %Fun 4
ICfuns{end+1}=@(x,y) moll(x,y,0,0,0.31,0.31);
    %Strain 4 patch
S=  [-0.4515, 0.4968, -0.9643, 0.3418];
dx= [-0.6988, 1.4363, -0.1722, -1.5009];
dy= [-1.7756, -1.4566, 0.4175, -0.0937];
p=  [0.6768, 0.3294, 0.5807, 0.2504];
for m=1:4
    ICfuns{end+1}=@(x,y) S(m)*exp(-((x-dx(m)).^2+(y-dy(m)).^2)./p(m)^2);
end

%Iterate over each of the IC funs
for IC=5:8%numel(ICfuns)
    w=w+ICfuns{IC}(wxm,wym);
end

%Solver--------------------------------------------------------------------
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
srcx=wxm(Nnumy(:,1,Estreamy(:,end))); %Global source location
srcy=wym(Nnumy(:,1,Estreamy(:,end)));
%Calculate global kernel for boundary velocity points
gkernel_xB= permute(bsxfun(@minus,rv_xB(1,2,:,:),srcy(:))./(sum(bsxfun(@minus,rv_xB,[srcx(:),srcy(:)]).^2,2)+2*del^2).^(3/2),[1 4 3 2]);
gkernel_yB= permute(bsxfun(@minus,srcx(:),rv_yB(1,1,:,:))./(sum(bsxfun(@minus,rv_yB,[srcx(:),srcy(:)]).^2,2)+2*del^2).^(3/2),[1 4 3 2]);
gkernel_x= squeeze(bsxfun(@minus,rv_x(1,2,:,:),srcy(:))./(sum(bsxfun(@minus,rv_x,[srcx(:),srcy(:)]).^2,2)+2*del^2).^(3/2));
gkernel_y= squeeze(bsxfun(@minus,srcx(:),rv_y(1,1,:,:))./(sum(bsxfun(@minus,rv_y,[srcx(:),srcy(:)]).^2,2)+2*del^2).^(3/2));

%Outer product of vorticity quadrature weights for pre-multiplication,
%including Jacobian
QwPre=(delX/2)^2*reshape(Qw'*Qw,1,[]);
k2=   zeros(size(wx));          %LSERK stage state
mask=0;
itt=0;
w_tot_elem=0;
lap=0;
lapper=1;
zmax=1.5*max(max(w));
zmin=1.5*min(min(w));
tic
for t=0:delt:endtime
    if mod(t,skip*delt)<delt
        itt=itt+1;
        tt(itt)=t;
        wxt(:,:,:,itt)=wx;
        surf(wxm,wym,reshape(wx,Np*K(1),Np*K(2))')
        axis([B,zmin,zmax])
        %Residual calc, used to calc the L^2 norm
        %GA1*exp(-((mod((wxm+1)/2-t*cx/2,1)*2-1).^2/Ga1+(mod((wym+1)/2-t*cy/2,1)*2-1).^2/Gb1))
        R=w-reshape(wx,Np*K(1),Np*K(2))';
        text(B(1),B(4),zmax*1.5*1.2,['Time: ',num2str(t),char(10),...
            'L^2 norm: ',num2str(sqrt(sum(sum(norm_h.*R.^2)))),char(10),...
            'Mask: ',num2str(length(mask)),char(10),...
            'Done in: ',num2str((endtime-t)*toc/t),char(10),...
            '\omega_{tot}: ',num2str(sum(w_tot_elem))]);
        text(B(1),B(4),-.5,num2str(lap));
        pause(0.0001)
    end
%     if mod(t,12)<eps(1)
%             lap(lapper)=sqrt(sum(sum(norm_h.*R.^2)))
%             lapper=lapper+1;
%     end

    %---Velocity eval of current timestep's vorticity config-----------
        v_xB(:)=0; v_yB(:)=0; v_xI(:)=0; v_yI(:)=0;
        w_elem=reshape(permute(reshape(wy,Np,K(2),Np,K(1)),[1 3 2 4]),1,Np^2,K(2)*K(1)); %Reshaped to col-wise element chunks
        w_tot_elem=abs(permute(mtimesx(w_elem,QwPre'),[3 1 2])); %Sum of vorticity in each elem
        mask=find(w_tot_elem>w_thresh); %Find "important" elements
        Nmask=find(not(w_tot_elem>w_thresh));
        w_elemPre=bsxfun(@times,QwPre,w_elem(:,:,mask)); %Pre-multiply by quad weights for speed
        
        msx=sort(reshape(Estreamx(:,mask),[],1));
        msy=sort(reshape(Estreamy(:,mask),[],1));
        Nmsx=sort(reshape(Estreamx(:,Nmask),[],1));
        Nmsy=sort(reshape(Estreamy(:,Nmask),[],1));
        for it=1:length(mask)
            w_source=w_elemPre(:,:,it);
            source= mask(it);
            
            %Form local source kernel from local stencil applied to global source kernel
            kernel_xB= gkernel_xB(:, [1:Np*K(2)] +Np*(K(2)-Enumy(source)), [1:K(1)+1] +(K(1)-Enumx(source)) );
            kernel_yB= gkernel_yB(:, [1:Np*K(1)] +Np*(K(1)-Enumx(source)), [1:K(2)+1] +(K(2)-Enumy(source)) );
            v_xB= v_xB + permute(mtimesx(w_source,kernel_xB),[2 3 1]);
            v_yB= v_yB + permute(mtimesx(w_source,kernel_yB),[3 2 1]);
            
            if Mp>2
                kernel_x= reshape(gkernel_x(:, [1:K(1)*(Mp-2)] +(Mp-2)*(K(1)-Enumx(source)), [1:Np*K(2)] +Np*(K(2)-Enumy(source)) ),Np^2,Mp-2,[]);
                kernel_y= reshape(gkernel_y(:, [1:K(2)*(Mp-2)] +(Mp-2)*(K(2)-Enumy(source)), [1:Np*K(1)] +Np*(K(1)-Enumx(source)) ),Np^2,Mp-2,[]);
                v_xI= v_xI + mtimesx(w_source,kernel_x);
                v_yI= v_yI + mtimesx(w_source,kernel_y);
            end
        end
        switch DGmask
            case 'reduced'
            %Assemble elementwise velocities in mask
            v_xE=[v_xB(EBl(msx)),v_xI(:,:,msx),v_xB(EBr(msx))];
            v_yE=[v_yB(EBb(msy)),v_yI(:,:,msy),v_yB(EBt(msy))];
            %Assemble boundary velocities not in mask
            v_xBQ=[v_xB(EBl(Nmsx)),v_xB(EBr(Nmsx))];
            v_yBQ=[v_yB(EBb(Nmsy)),v_yB(EBt(Nmsy))];
            case 'full'
            v_xE=[v_xB(EBl),v_xI,v_xB(EBr)];
            v_yE=[v_yB(EBb),v_yI,v_yB(EBt)];
        end
        %---Velocity eval ends---------------------------------------------
        
    for i=1:nS
        St= t+RKc(i)*delt;              %Unused currently, St is the stage time if needed
        
        
        
        %---Advection------------------------------------------------------
        w_lx= mtimesx(Ll',wx);          %Left interpolated vorticity
        w_rx= mtimesx(Lr',wx);          %Right interpolated vorticity
        w_bx= mtimesx(Ll',wy);          %Bottom interpolated vorticity
        w_tx= mtimesx(Lr',wy);          %Top interpolated vorticity
        
        %Boundary fluxes
        if BCtype== 'NoInflow'
            v_xBC= v_xB; v_xBC(:,1)= min(v_xBC(:,1),0); v_xBC(:,end)= max(v_xBC(:,end),0);
            v_yBC= v_yB; v_yBC(1,:)= min(v_yBC(1,:),0); v_yBC(end,:)= max(v_yBC(end,:),0);
        end
        fl= abs( v_xBC(EBl) ).*( w_rx(x_km1).*(sign(v_xB(EBl))+alpha) + w_lx.*(sign(v_xB(EBl))-alpha) );
        fr= abs( v_xBC(EBr) ).*( w_rx.*(sign(v_xB(EBr))+alpha) + w_lx(x_kp1).*(sign(v_xB(EBr))-alpha) );
        fb= abs( v_yBC(EBb) ).*( w_tx(y_km1).*(sign(v_yB(EBb))+alpha) + w_bx.*(sign(v_yB(EBb))-alpha) );
        ft= abs( v_yBC(EBt) ).*( w_tx.*(sign(v_yB(EBt))+alpha) + w_bx(y_kp1).*(sign(v_yB(EBt))-alpha) );
        
        %Nodal total surface flux
        SurfFlux_x=bsxfun(@times,fr,LrM)-bsxfun(@times,fl,LlM);
        SurfFlux_y=bsxfun(@times,ft,LrM)-bsxfun(@times,fb,LlM);
        %Nodal stiffness eval
        switch DGmask
            case 'reduced'
            Stiff_x(:,:,msx,:)= mtimesx(v_xE,mtimesx(QwSM,wx(:,:,msx)));
            Stiff_x(:,:,Nmsx,:)= mtimesx(v_xBQ,mtimesx(QwSMlow,wx(:,:,Nmsx)));
            Stiff_y(:,:,msy,:)= mtimesx(v_yE,mtimesx(QwSM,wy(:,:,msy)));
            Stiff_y(:,:,Nmsy,:)= mtimesx(v_yBQ,mtimesx(QwSMlow,wy(:,:,Nmsy)));
            case 'full'
            Stiff_x= mtimesx(v_xE,mtimesx(QwSM,wx));
            Stiff_y= mtimesx(v_yE,mtimesx(QwSM,wy));
        end

        wx_dt= permute(Stiff_x-SurfFlux_x,[4 1 3 2]); %Reshape to match wx
        wy_dt= reshape(reshape(Stiff_y-SurfFlux_y,K(2),[])',Np,1,[]); %Reshape to match wx
        
        k2= RKa(i)*k2 + delt*(wx_dt+wy_dt);
        wx= wx+RKb(i)*k2;
        wy= reshape(reshape(wx,K(1)*Np,[])',Np,1,[]); %Reshape wx to match global node ordering
    end
end
a=sum(w_tot_elem);
wy= reshape(w,Np,1,[]);
w_elem=reshape(permute(reshape(wy,Np,K(2),Np,K(1)),[1 3 2 4]),1,Np^2,K(2)*K(1));
w_tot_elem=abs(permute(mtimesx(w_elem,QwPre'),[3 1 2]));
b=sum(w_tot_elem);
a-b