%--Josh Bevan 2014
%--1D Scalar Conservation Eqn solution using Discontinuous Galerkin
clear all
close all
tau=2*pi();

%--Here we attempt to solve the 1D scalar conservation eqn of the form:
%\frac{\partial u)(\partial t} + \frac{\partial f(u))(\partial x} = 0
%where f(u) is some flux function describing the "flow" of a conserved
%quantity. In this simplified case f(u) = u giving us a linear PDE
%We will use periodic boundary conditions to examine numerical dissipation
%effects

%--Let u(x,t) be the exact solution on the domain 0<=x<=1
%--Let u(.,0) = u0 = sin(2pi * x)

%--Let Vh be the finite vector space of linear polynomials of order N.
%--Let v be the approximate numerical solution consisting of a linear 
%combination of basis functions (\psiN for the Nth basis) in Vh with
%scalar coefficients BasisWeightN (i.e. \sum BasisWeightN \psiN)
%--Let \phiN be the Nth test function in the same vector space as the basis
%functions (Vh)
N=5;    

%--Discretize the domain into K elements with K+1 nodes, we use a constant
%spacing, but it could be arbitrary
K=16;
xNode=0:1/K:1;
deltax = diff(xNode);
%All the element boundaries, note there are repeats since nodes have a coincident
%brother from an adjacent element, except at the domain boundaries
elemBC = reshape(sort([xNode,xNode(2:end-1)]),2,K)';
u= sin(tau.*elemBC);

%--According to Cockburn,Shu 2001 Eqn 2.2 let uh(.,0) be computed by
%\int v \phi = \int u0 \phi for each element (xk-1/2 < x < xk+1/2)
%--We can explicitly define a formula for the value of the RHS 
%ExactRHSN = \int u0 \phiN
for k=1:length(elemBC)
    [x,map,w]=gauss(N,elemBC(k,1),elemBC(k,2));

    for m=0:(N+1)/2
            l = legendre(m,x);
            l = l(1,:);
            BasisWeights(m+1) = (m+.5)*sum(w.*sin(tau*map)'.*l);
    end
end