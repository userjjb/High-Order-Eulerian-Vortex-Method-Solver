function [W,Ll,Lr] = StiffnessQuadModWeights(Qx,Qx2)
%Returns the M,N,N matrix of modified quadrature weights associated with
%the N quadrature nodes Qx and M nodes Qx2 for each weighting basis m<=M.
%Qx should be the conserved quantity and Qx2 the flux velocity. This 
%quadrature is used to produce a stiffness matrix that is the integral of
%the product of two interpolations times the spatial deriavitve of the
%weighting function.
%This is distinct from the integral of the interpolation of the product of
%two functions times the derivative of the weighting function, and more
%accurate. It also permits the use of mixed quadrature, i.e. different
%order quadrature for each interpolation as well as differnet nodes (GL
%mixed with LGL for instance).
N=length(Qx);
M=length(Qx2);

nn=elim(Qx(1:N)',Qx(1:N)',[1 3 2]); %All points except the basis point for each basis
Lag1= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,nn(nv,:,:)),bsxfun(@minus,Qx(nv),nn(nv,:,:))),3); %Associated Lagrange basis function
nn2=elim(Qx2(1:M)',Qx2(1:M)',[1 3 2]);
Lag2= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,nn2(nv,:,:)),bsxfun(@minus,Qx2(nv),nn2(nv,:,:))),3);
dLag1= @(x,nv) Lag1(x,nv).*sum(1./bsxfun(@minus,x,nn(nv,:,:)),3);

W=zeros(M,N,1,N);
for j=1:N
    for i=1:M
        for n=1:N
            trip=@(x) Lag1(x,j).*Lag2(x,i).*dLag1(x,n);
            W(i,j,1,n)=integral(trip,-1,1,'RelTol',1e-14,'AbsTol',1e-15);
        end
    end
end

Ll=Lag1(-1,1:N);
Lr=Lag1(1,1:N);
end