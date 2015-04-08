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
N= 10;                               %Local vorticity poly order
M= N;                               %Local velocity poly order
%Global domain initialization (parameters)---------------------------------
B= [-1 1 -1 1];                     %left, right, bottom, top
K= [10 10];                         %Num elements along x,y
Ex= linspace(B(1),B(2),K(1)+1);     %Elem edges left-right
Ey= linspace(B(3),B(4),K(2)+1);     %Elem edges bottom-top
%Global domain initialization (calculated)---------------------------------
Np= N+1;                            %Number of interpolation points to achieve N
Mp= M+1;                            %Number of interpolation points to achieve M
delX= (B(2)-B(1))/K(1);             %Element size
    %For simplicty we require square elements
    assert(abs(delX-(B(4)-B(3))/K(2))<eps(delX),'Elements must be square')
    
[Qx, Qw]= GLquad(Np);               %Gauss-Legendre quadrature nodes/weights
[Qx2, Qw2]= GLquad(Mp);             %
[QwS, Ll, Lr]= StiffnessQuadModWeights(Qx,Qx2); %Modified quadrature weights for the stiffness matrix calculation, and vorticity interpolation evals
QwSM=bsxfun(@times,QwS,2./(delX*permute(Qw,[1 4 3 2]))); %Distributed Mass matrix and Jacobian

w= zeros(Np*K(2),Np*K(1));          %Global vorticity at element interp points
v_x= ones(Mp*K(2),Mp*K(1));         %Global velocity_x at element interp points
v_y= zeros(Mp*K(2),Mp*K(1));        %Global velocity_y at element interp points
v_xB= ones(Mp*K(2),K(1)+1);         %Global velocity_x at elem x_boundaries
v_yB= zeros(K(2)+1,Mp*K(1));        %Global velocity_y at elem y_boundaries

%Element global numbering as well as element to left or right for periodic
%BCs For instance,pass x_km1 as an index to get the element left of current
eord_x= reshape(1:K(1)*K(2)*Np,K(1),K(2)*Np)';
x_km1= reshape(circshift(eord_x,[0 1])',[],1);
x_kp1= reshape(circshift(eord_x,[0 -1])',[],1);
eord_y= reshape(1:K(1)*K(2)*Np,K(2),K(1)*Np);
y_km1= reshape(circshift(eord_y,[1 0])',[],1);
y_kp1= reshape(circshift(eord_y,[-1 0])',[],1);
%Element bound numbering
EBl= reshape(bsxfun(@plus,1:K(1),(K(1)+1)*[0:Mp*K(2)-1]')',1,1,[]);
EBr= reshape(bsxfun(@plus,2:K(1)+1,(K(1)+1)*[0:Mp*K(2)-1]')',1,1,[]);
EBt= reshape(bsxfun(@plus,1:K(2),(K(2)+1)*[0:Mp*K(1)-1]'),1,1,[]);
EBb= reshape(bsxfun(@plus,2:K(2),(K(2)+1)*[0:Mp*K(1)-1]'),1,1,[]);

%Interp/quad node locations
x_w= reshape(bsxfun(@plus,repmat((Qx+1)*(delX/2),1,K(1)),Ex(1:end-1)),[],1);
y_w= reshape(bsxfun(@plus,repmat((Qx+1)*(delX/2),1,K(2)),Ey(1:end-1)),[],1);
%Elementwise velocity interp/quad node locations
x_v= reshape(bsxfun(@plus,repmat((Qx2+1)*(delX/2),1,K(1)),Ex(1:end-1)),[],1);
y_v= reshape(bsxfun(@plus,repmat((Qx2+1)*(delX/2),1,K(2)),Ey(1:end-1)),[],1);
%Note that the actual velocity node locations for x_streams are (x_v,y_w) and
%for y_streams are (y_v,x_w). The velocity grid is NOT a tensor product, but
%rather velocity "streams" are coincident with vorticity streams, to reuse
%evaluations for both x and y streams we set Mp=Np so we get a tensor grid
%as a byproduct, but it isn't required
%Additionally, the boundary velocity locations occur at (Ex,y_w) and
%(x_w,Ey) analogous to the interpolation points
[wxm, wym]= meshgrid(x_w,y_w);

%Intial conditions specification-------------------------------------------
ICfuns={}; %Create a cell list of functions that define the ICs
    %Gaussian 1
Ga1=    0.2;
Gb1=    0.05;
Gdx1=   0; 
Gdy1=   0;
GA1=    1;
ICfuns{end+1}=@(x,y) GA1*exp(-((x-Gdx1).^2/Ga1+(y-Gdy1).^2/Gb1));%Center is at (dx,dy)

%Iterate over each of the IC funs
for IC=1:numel(ICfuns)
    w=w+ICfuns{IC}(wxm,wym);
end
%clearvars wxm wym %No longer need the meshgrid points

%Solver--------------------------------------------------------------------
%Semi-discrete-------------------------------------------------------------
%Setup

%Weighting basis with lumped outside terms to save on calculation for
%numerical total nodal surface flux term: \hat{f}_R L_j(x_R)-\hat{f}_L L_j(x_L)
%Distributed: inv Mass matrix (1/QwG), inv Jacobian (2/delX), and 1/2 coeff
%from num flux eval
LlM=  permute( Ll.*(1./(delX*Qw')) , [2 3 4 1]); %Left elem boundary eval
LrM=  permute( Lr.*(1./(delX*Qw')) , [2 3 4 1]); %Right elem boundary eval

wx=   reshape(w',Np,1,[]);      %Reshape vorticity for mtimesx bsx
v_x=  reshape(v_x',1,Mp,[]);    %Reshape velocity_x for mtimesx bsx

%Loop starts here
delt=0.001;
skip=0.05;
for t=0:delt:3
    w_lx= mtimesx(Ll',wx);          %Left interpolated vorticity
    w_rx= mtimesx(Lr',wx);          %Right interpolated vorticity
    %Boundary fluxes
    fl= abs( v_x(EBl) ).*( w_rx(x_km1).*(sign(v_x(EBl))+alpha) + w_lx.*(sign(v_x(EBl))-alpha) );
    fr= abs( v_x(EBr) ).*( w_rx.*(sign(v_x(EBr))+alpha) + w_lx(x_kp1).*(sign(v_x(EBr))-alpha) );

    SurfFlux_x=bsxfun(@times,fr,LrM)-bsxfun(@times,fl,LlM); %Nodal total surface flux
    Stiff_x= mtimesx(v_x,mtimesx(QwSM,wx)); %Nodal stiffness eval

    wx_dt= permute(Stiff_x-SurfFlux_x,[4 1 3 2]);
    wx= wx+wx_dt*delt;
    if mod(t,skip)<delt
        surf(wxm,wym,reshape(wx,Np*K(1),Np*K(2))')
        pause(0.001)
    end
end