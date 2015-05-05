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
v_xE=zeros(1,Mp,K(1)*Np*K(2));
v_yE=zeros(1,Mp,K(2)*Np*K(1));
v_xB=zeros(Np*K(2),K(1)+1);
v_yB=zeros(K(2)+1,Np*K(1));
v_xBF=zeros(1,2,Np*K(2)*K(1));
v_yBF=zeros(1,2,Np*K(1)*K(2));

%---Node numbering---------------------------------------------------------
%Element numbering
[Enumx,Enumy]= meshgrid(1:K(1),1:K(2));
Enum= reshape(1:K(1)*K(2),K(2),K(1));
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
Ex= linspace(B(1),B(2),K(1)+1);     %Elem edges left-right
Ey= linspace(B(3),B(4),K(2)+1);     %Elem edges bottom-top
%Interp/quad node locations
x_w= reshape(bsxfun(@plus,repmat((Qx+1)*(delX/2),1,K(1)),Ex(1:end-1)),[],1);
y_w= reshape(bsxfun(@plus,repmat((Qx+1)*(delX/2),1,K(2)),Ey(1:end-1)),[],1);
%Boundary velocity interp/quad node locations, reflected about both axes
if 1
    [t2,t1]= meshgrid([y_w;(B(4)-B(3))+y_w],[Ex, (B(2)-B(1))+Ex(2:end)]);
    rv_xB= reshape([t1(:),t2(:)]',1,2,2*K(1)+1,[]);
    [t1,t2]= meshgrid([x_w;(B(2)-B(1))+x_w],[Ey,(B(4)-B(3))+Ey(2:end)]);
    rv_yB= reshape([t1(:),t2(:)]',1,2,2*K(2)+1,[]);
else
    [t2,t1]= meshgrid(y_w,Ex);
    rv_xB= reshape([t1(:),t2(:)]',1,2,2*K(1)+1,[]);
    [t1,t2]= meshgrid(x_w,Ey);
    rv_yB= reshape([t1(:),t2(:)]',1,2,2*K(2)+1,[]);
end

%Internal element interp/quad node locations
x_v=reshape(bsxfun(@plus,(Qx2(2:end-1)+1)*(delX/2),Ex(1:end-1)),1,[]);
y_v=reshape(bsxfun(@plus,(Qx2(2:end-1)+1)*(delX/2),Ey(1:end-1)),1,[]);
[t2,t1]= meshgrid(y_w,x_v);
rv_x= reshape([t1(:),t2(:)]',1,2,Mp-2,[]);
[t1,t2]= meshgrid(x_w,y_v);
rv_y= reshape([t1(:),t2(:)]',1,2,Mp-2,[]);
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

%Tidy up un-needed vars
clearvars h_x h_y t1 t2 