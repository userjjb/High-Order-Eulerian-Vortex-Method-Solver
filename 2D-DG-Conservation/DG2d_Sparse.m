close all
clear all
clc

%Global domain initialization (parameters)---------------------------------
B= [-1 1 -1 1];                     %left, right, bottom, top
N= 4;                               %Local vorticity poly order
M= 4;                               %Local velocity poly order
K= [10 10];                         %Num elements along x,y
Ex= linspace(B(1),B(2),K(1)+1);     %Elem edges left-right
Ey= linspace(B(3),B(4),K(2)+1);     %Elem edges bottom-top
%Global domain initialization (calculated)---------------------------------
Np= N+1;                            %Number of interpolation points to achieve N
Mp= M+1;                            %Number of interpolation points to achieve M
delX= (B(2)-B(1))/K(1);             %Element size
    %For simplicty we require square elements
    assert(abs(delX-(B(4)-B(3))/K(2))<eps(delX),'Elements must be square')
    
[QxG, QwG]= GLquad(Np);             %Gauss-Legendre quadrature nodes/weights
[QxL, QwL]= LGLquad(Mp);            %Lobatto-Gauss-Legendre quadrature nodes/weights
[QwS, Ll, Lr]= StiffnessQuadModWeights(QxG,QxL); %Modified quadrature weights for the stiffness matrix calculation, and vorticity interpolation evals
QwSM=bsxfun(@times,QwS,2./(delX*permute(QwG,[1 4 3 2]))); %Distributed Mass matrix and Jacobian

w= zeros(Np*K(2),Np*K(1));          %Global vorticity at interp points
v_x= ones(Mp*K(2),Mp*K(1));         %Global velocity_x at interp points
v_y= zeros(Mp*K(2),Mp*K(1));        %Global velocity_y at interp points
ord_x= reshape(1:Np^2*K(2)*K(1),Np*K(1),Np*K(2))'; %Global node numbering order along stream_x
ord_y= reshape(1:Np^2*K(2)*K(1),Np*K(2),Np*K(1)); %Global node numbering order along stream_y

%Interp/quad node locations
wx= reshape(bsxfun(@plus,repmat((QxG+1)*(delX/2),1,K(1)),Ex(1:end-1)),[],1);
wy= reshape(bsxfun(@plus,repmat((QxG+1)*(delX/2),1,K(2)),Ey(1:end-1)),[],1);
%Global velocity interp/quad node locations
vxG= [reshape(bsxfun(@plus,repmat((QxL(1:end-1)+1)*(delX/2),1,K(1)),Ex(1:end-1)),[],1) ; B(2)];
vyG= [reshape(bsxfun(@plus,repmat((QxL(1:end-1)+1)*(delX/2),1,K(2)),Ey(1:end-1)),[],1) ; B(4)];
%Elementwise velocity interp/quad node locations
vxE= reshape(bsxfun(@plus,repmat((QxL+1)*(delX/2),1,K(1)),Ex(1:end-1)),[],1);
vyE= reshape(bsxfun(@plus,repmat((QxL+1)*(delX/2),1,K(2)),Ey(1:end-1)),[],1);
[wxm, wym]=meshgrid(wx,wy);

%Intial conditions specification-------------------------------------------
ICfuns={}; %Create a cell list of functions that define the ICs
    %Gaussian 1
Ga1=0.2;
Gb1=0.05;
Gdx1=0; 
Gdy1=0;
GA1=1;
ICfuns{end+1}=@(x,y) GA1*exp(-((x-Gdx1).^2/Ga1+(y-Gdy1).^2/Gb1));%Center is at (dx,dy)

%Iterate over each of the IC funs
for IC=1:numel(ICfuns)
    w=w+ICfuns{IC}(wxm,wym);
end
clearvars wxm wym %No longer need the meshgrid points

%Solver--------------------------------------------------------------------
%Semi-discrete-------------------------------------------------------------
%Setup
LlM= Ll.*(2./(delX*QwG'));                  %Distributed Mass matrix and Jacobian
LrM= Lr.*(2./(delX*QwG'));                  %Distributed Mass matrix and Jacobian

w= permute(reshape(w',Np,[]),[1 3 2]);      %Reshape vorticity for mtimesx bsx
v_x= permute(reshape(v_x',Mp,[]),[3 1 2]);  %Reshape velocity_x for mtimesx bsx

w_km1=w(reshape(circshift(ord_x,[0 Np])',Np,[]));%FIX ME NOT RIGHT COMAPARE TO PROPER WRAP for w_k-1

Stiff_x= permute(mtimesx(v_x,mtimesx(QwSM,w)),[3 4 1 2]);
NumFlux_x=
Boundary_x= bsxfun(@times,,LrM)-bsxfun(@times,,LlM);













