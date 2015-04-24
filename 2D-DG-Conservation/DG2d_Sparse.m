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
N= 4;                               %Local vorticity poly order
M= 4;                               %Local velocity poly order
[RKa,RKb,RKc,nS]= LSRKcoeffs('NRK14C');
w_thresh=0;
%---Global domain initialization (parameters)------------------------------
B= [-1 1 -1 1];                     %left, right, bottom, top
K= [16 16];                         %Num elements along x,y
Ex= linspace(B(1),B(2),K(1)+1);     %Elem edges left-right
Ey= linspace(B(3),B(4),K(2)+1);     %Elem edges bottom-top

c=1;
theta = pi/4;
cx = c*cos(theta);
cy = c*sin(theta);
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
v_xI=cx*ones(Mp-2,1,K(1)*Np*K(2));
v_yI=cy*ones(Mp-2,1,K(2)*Np*K(1));
v_xB=cx*ones(Np*K(2),K(1)+1);
v_yB=cy*ones(K(2)+1,Np*K(2));

%---Node numbering---------------------------------------------------------
%Element global numbering as well as element to left or right for periodic
%BCs For instance,pass x_km1 as an index to get the element left of current
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
EBy=reshape(1:Np*K(1)*(K(2)+1),(K(1)+1),Np*K(1));
EBb= reshape(EBy(1:end-1,:),1,1,[]);
EBt= reshape(EBy(2:end,:),1,1,[]);
%Stream/element associativity (stream in:,element) col-wise
Estreamx= reshape(Eord_x,Np,[]);
Estreamy= reshape(permute(reshape(Eord_y',Np,K(1),[]),[1 3 2]),Mp,[]);
%Node/element associativity

%---Node coordinates-------------------------------------------------------
%Interp/quad node locations
x_w= reshape(bsxfun(@plus,repmat((Qx+1)*(delX/2),1,K(1)),Ex(1:end-1)),[],1);
y_w= reshape(bsxfun(@plus,repmat((Qx+1)*(delX/2),1,K(2)),Ey(1:end-1)),[],1);
%Elementwise velocity interp/quad node locations
x_v= reshape(bsxfun(@plus,repmat((Qx2+1)*(delX/2),1,K(1)),Ex(1:end-1)),[],1);
y_v= reshape(bsxfun(@plus,repmat((Qx2+1)*(delX/2),1,K(2)),Ey(1:end-1)),[],1);
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
moll=@(x,y,dx,dy,a,b) heaviside(1-(((x-dx)/a).^2+((y-dy)/b).^2)).*exp(1+(-1./(1-(((x-dx)/a).^2+((y-dy)/b).^2).^4)));
%Intial conditions specification-------------------------------------------
ICfuns={}; %Create a cell list of functions that define the ICs
    %Gaussian 1
Ga1=    0.05;
Gb1=    0.05;
Gdx1=   0; 
Gdy1=   0;
GA1=    1;
ICfuns{end+1}=@(x,y) GA1*exp(-(((x)-Gdx1).^2/Ga1+(y-Gdy1).^2/Gb1));%Center is at (dx,dy)
    %Gaussian 2
Ga2=    0.02;
Gb2=    0.05;
Gdx2=   0.5; 
Gdy2=   -0.7;
GA2=    .4;
ICfuns{end+1}=@(x,y) GA2*exp(-(((x)-Gdx2).^2/Ga2+(y-Gdy2).^2/Gb2));%Center is at (dx,dy)

%Iterate over each of the IC funs
for IC=1:1%numel(ICfuns)
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

%Intialize scalar kernel values, a specific stream's kernel should be a
%matrix of pointwise values: each column corresponds to the matching point
%in the stream, each element in the column matches the scalar value of the
%kernel for that point as one moves col-wise within the source
%element to match w_elem's ordering.
gkernelxB= ;
gkernelyB= ;
gkernelx= ;
gkernely= ;

oner=ones(Np^2,1);              %For fast element sums
QwPre=reshape(Qw'*Qw,1,[]);     %Outer product of vorticity quadrature weights for pre-multiplication
delt= 0.04;
skip= 0.1;
k2=   zeros(size(wx));          %LSERK stage state
for t=0:delt:10
    if mod(t,skip)<delt
        surf(wxm,wym,reshape(wx,Np*K(1),Np*K(2))')
        axis([B,0,1.1])
        text(-1,1,1.2,strcat('Time: ',num2str(t)));
        %Residual calc, used to calc the L^2 norm
        R=GA1*exp(-((mod((wxm+1)/2-t*cx/2,1)*2-1).^2/Ga1+(mod((wym+1)/2-t*cy/2,1)*2-1).^2/Gb1))-reshape(wx,Np*K(1),Np*K(2))';
        text(-1,1,1.3,strcat('L^2 norm: ',num2str(sqrt(sum(sum(norm_h.*R.^2))))));
        pause(0.0001)
    end
    for i=1:nS
        St= t+RKc(i)*delt;              %Unused currently, St is the stage time if needed
        %---Velocity eval of current timestep's vorticity config-----------
        v_xB(:)=0; v_yB(:)=0; v_xI(:)=0; v_yI(:)=0;
        w_elem=reshape(permute(reshape(wy,Np,K(2),Np,K(1)),[1 3 2 4]),1,Np^2,K(2)*K(1)); %Reshaped to col-wise element chunks
        w_tot_elem=abs(permute(mtimesx(w_elem,oner),[3 1 2])); %Sum of vorticity in each elem
        mask_elem=find(w_tot_elem>w_thresh); %Find "important" elements
        w_elemPre=bsxfun(@times,QwPre,w_elem(:,:,mask_elem)); %Pre-multiply by quad weights for speed
        
        mask_streamx=reshape(Estreamx(:,mask_elem),1,1,[]);
        mask_streamy=reshape(Estreamy(:,mask_elem),1,1,[]);
        for source=1:length(mask_elem)
            %Element mask transformed to global coords
            mask_tfB= ;
            mask_tf= ;
            w_source=w_elemPre(:,:,source);
            %Form the (Np^2,Mp-2 OR 2,Np*length(mask_elem)) element kernel
            %size comes from(num w_interp,velocity eval,streams in total)
            kernel_xB=gkernel_xB(:,:,mask_tfB);
            kernel_yB=gkernel_yB(:,:,mask_tfB);
            kernel_x=gkernel_x(:,:,mask_tf);
            kernel_y=gkernel_y(:,:,mask_tf);
            v_xB= v_xB + mtimesx(w_source,kernelxB);
            v_yB= v_yB + mtimesx(w_source,kernelyB);
            v_xI(mask_streamx)= v_xI(mask_streamx) + mtimesx(w_source,kernelx);
            v_yI(mask_streamy)= v_yI(mask_streamy) + mtimesx(w_source,kernely);
        end
        v_xE=[v_xB(EBl),v_xI,v_xB(EBr)]; %Calc elementwise velocities
        v_yE=[v_yB(EBb),v_yI,v_yB(EBt)];
        
        %---Advection------------------------------------------------------
        w_lx= mtimesx(Ll',wx);          %Left interpolated vorticity
        w_rx= mtimesx(Lr',wx);          %Right interpolated vorticity
        w_bx= mtimesx(Ll',wy);          %Bottom interpolated vorticity
        w_tx= mtimesx(Lr',wy);          %Top interpolated vorticity
        %Boundary fluxes
        fl= abs( v_xB(EBl) ).*( w_rx(x_km1).*(sign(v_xB(EBl))+alpha) + w_lx.*(sign(v_xB(EBl))-alpha) );
        fr= abs( v_xB(EBr) ).*( w_rx.*(sign(v_xB(EBr))+alpha) + w_lx(x_kp1).*(sign(v_xB(EBr))-alpha) );
        fb= abs( v_yB(EBb) ).*( w_tx(y_km1).*(sign(v_yB(EBb))+alpha) + w_bx.*(sign(v_yB(EBb))-alpha) );
        ft= abs( v_yB(EBt) ).*( w_tx.*(sign(v_yB(EBt))+alpha) + w_bx(y_kp1).*(sign(v_yB(EBt))-alpha) );

        SurfFlux_x=bsxfun(@times,fr,LrM)-bsxfun(@times,fl,LlM); %Nodal total surface flux
        Stiff_x= mtimesx(v_xE,mtimesx(QwSM,wx)); %Nodal stiffness eval
        %Stiff_x()= mtimesx(v_xB,mtimesx(QwSMlow,wx)); %Nodal stiffness eval
        SurfFlux_y=bsxfun(@times,ft,LrM)-bsxfun(@times,fb,LlM);
        Stiff_y= mtimesx(v_yE,mtimesx(QwSM,wy));
        %Stiff_y()= mtimesx(v_yB,mtimesx(QwSMlow,wy));

        wx_dt= permute(Stiff_x-SurfFlux_x,[4 1 3 2]); %Reshape to match wx
        wy_dt= reshape(reshape(Stiff_y-SurfFlux_y,K(2),[])',Np,1,[]); %Reshape to match wx
        
        k2= RKa(i)*k2 + delt*(wx_dt+wy_dt);
        wx= wx+RKb(i)*k2;
        wy= reshape(reshape(wx,K(1)*Np,[])',Np,1,[]); %Reshape wx to match global node ordering
    end
end