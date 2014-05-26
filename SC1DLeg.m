%--Josh Bevan 2014
%--1D Scalar Conservation Eqn solution using Discontinuous Galerkin
close all
clear all

tau=2*pi();
N=3;
K=32;

%--Here we attempt to solve the 1D scalar conservation eqn of the form:
%\frac{\partial u)(\partial t} + \frac{\partial f(u))(\partial x} = 0
%where f(u) is some flux function describing the "flow" of a conserved
%quantity. In this simplified case f(u) = u giving us a linear PDE
%We will use periodic boundary conditions to examine numerical dissipation
%effects

%--Let u(x,t) be the exact solution on the domain 0<=x<=1
%--Let u(.,0) = u0 = sin(2pi * x)

%--Let Vh be the finite vector space of the shifted Legendre polynomials
%up to order N.
%--Let v be the approximate numerical solution consisting of a linear 
%combination of basis functions (\psiN for the Nth basis) in Vh with
%scalar coefficients BasisWeight_n (i.e. \sum_{n=0}^N BasisWeight_n \psi_n)
%--Let \phi_n be the nth test function in the same vector space as the basis
%functions (Vh)   

%--Discretize the domain into K elements with K+1 nodes, we use a constant
%spacing, but it could be arbitrary
xNode=0:1/K:1;
deltax = diff(xNode);
%All the element boundaries, note there are repeats since nodes have a coincident
%brother from an adjacent element, except at the domain boundaries
elemBC = reshape(sort([xNode,xNode(2:end-1)]),2,K)';

%--According to Cockburn,Shu 2001 Eqn 2.2 let uh(.,0) be computed by
%\int v \phi = \int u0 \phi for each element (xk-1/2 < x < xk+1/2)
%--We will numerically compute the RHS integral using Gauss-Legendre
%quadrature which is accurate up to 2N-1

%Precalculate quadrature nodes and weights 
[Qx,Qw]=gauss(N+1);

%Precompute Legendre values for fixed order quadrature points
L = zeros(N+1);
for m=0:N
    temp = legendre(m,Qx); %legendre() actually calculates the associated
    L(:,m+1) = temp(1,:)'; %legendre polys, we only need m=0 (first row)
end

%--We can algebraically solve for the basis weights by analytically computing
%the LHS integral. The integrand is the inner product of v0 and the test
%function.
%--Since both are orthogonal in Vh the result is a constant times
%the Kronecker product. The RHS is the discrete Legendre transform
%resulting from the Gauss-Legendre quadrature.
%--The simplified result is \overset{m}{a} =
%(m+.5) \sum_{i=1}^N Qw_i sin(2\pi \tilde{x}(x) ) \overset{n}{\phi}(x_i)
BasisWeights = zeros(K,N+1);
map = zeros(K,N+1);
for k=1:K
    map(k,:) = (elemBC(k,2)-elemBC(k,1))*Qx/2 + (elemBC(k,1)+elemBC(k,2))/2;
    for m=0:N
        BasisWeights(k,m+1) = (m+.5)*sum(Qw.*sin(tau*map(k,:)).*L(:,m+1)');
    end
end

%--Now that we have u0 we can begin explicit time stepping (forward Euler)
%with the semi-discrete form of the PDE.
%--Because we have a linear flux function all the classic monotone flux
%schemes reduce to the simple upwind flux i.e. g(v-(x),g(v+(x))=v-(x)

%Self-references to the bases in the element. This is the sum of part of
%the upwind flux and the RHS integral in the original PDE. The RHS integral
%portion is a strictly lower triangular matrix
SelfStencil = 2*(toeplitz(mod(0:N,2),0:N<0))-ones(N+1);

%References to the bases in the upwind element, consists of the other part
%of the upwind flux
UpwindStencil = ones(N+1); UpwindStencil(2:2:N+1,:)=-1;

%LHSIntegral is a purely diagonal matrix that depends on element size and
%the order of each basis
LHSIntegral = reshape([1:2:2*N+1]'*(1./deltax),K*(N+1),1);

%Build a unit stencil to pull the global diagonal pattern from
C=SelfStencil; C(N+2:2*N+2,1:N+1)=UpwindStencil;
DiagStencil = zeros(N+1,3*N+2);
for d=-(2*N+1):N
    DiagStencil(1:length(diag(C,d)),d+2*N+2) = diag(C,d);
end

%Flip superdiagonal since spdiags truncates the start of any oversized
%superdiagonals
DiagStencil(:,2*N+2:3*N+2) = flipud(DiagStencil(:,2*N+2:3*N+2));

%Build sparse reference matrix by repeating unit stencil for each element
A = spdiags(repmat(DiagStencil,K,1),-(2*N+1):N,K*(N+1),K*(N+1));
A(1:N+1,K*(N+1)-N:K*(N+1))=UpwindStencil;
spy(A)

%Calculate new da/dt
%BasisWeights_dt=LHSIntegral.*(A*reshape(BasisWeights',K*(N+1),1));
BasisWeights = reshape(BasisWeights',K*(N+1),1);
    
%--Discretize in time and plot
deltaT= .0001;
saveT = 0.05;
endT = 10;

nsaveT = floor(saveT/deltaT);
nT = floor(endT/deltaT);
%saved = zeros(N+1,K,(nT/nsaveT)+1);
i=1;

saved(:,:,1)=reshape(BasisWeights,N+1,K);
for t= 1:1:nT
    if t/nsaveT==floor(t/nsaveT)    
        i= i+1;
        saved(:,:,i)=reshape(BasisWeights,N+1,K);
    end
    BasisWeights_dt=LHSIntegral.*(A*BasisWeights);
    BasisWeights = BasisWeights+(BasisWeights_dt.*deltaT);
end

Lobatto = zeros(N+3,N+1);
for m=0:N
    temp = legendre(m,[-1; Qx; 1]);
    Lobatto(:,m+1) = temp(1,:)';
end
j=0;
for i=1:length(saved)
    plot([elemBC(:,1) map elemBC(:,2)]',(Lobatto*saved(:,:,i)),'o-')
    axis([0 1 -1.5 1.5])
    text(1.02,-0.2,'RMS')
    text(1.02,-.3,num2str(rms(reshape(L*saved(:,:,i),K*(N+1),1))));
    text(1.02,1.1,'Time')
    text(1.02,1,num2str((i-1)*saveT));
    pause(.01)
end