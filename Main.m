clear all

%%Let u(x,t) be the exact solution on the domain 0<=x<=1
%%Let u(.,0) = u0 = sin(2pi * x)

%%Let vh be the finite vector space of linear polynomials
%%Let uh be the approximate numerical solution consisting of a linear 
%combination of basis functions (\thetaN for the Nth basis) in vh with
%scalar coefficients BasisWeightN (i.e. \sum BasisWeightN \thetaN)
%%Let \phiN be the Nth test function in the same vector space as the basis
%functions (vh)

%%Discretize the domain into K elements with K+1 nodes
tau=2*pi();

K=20;
xNode=0:1/K:1;
deltax = diff(xNode);
%All the nodes, note their are repeats since nodes have a coincident
%brother from an adjacent element, except at the boundaries
xh = sort([xNode,xNode(2:end-1)]);

xf = 0:1/100:1;
u= sin(tau.*xf);

%%According to Cockburn,Shu 2001 Eqn 2.2 let uh(.,0) be computed by
%\int uh \phi = \int u0 \phi for each element (xj-1/2 < x < xj=1/2)
%%We can explicitly define a formula for ExactRHSN = \int u0 \phiN the value of
%the RHS
ExactRHS(1,:)= (tau*deltax.*cos(tau.*xNode(1:end-1))+sin(tau.*xNode(1:end-1))-sin(tau.*xNode(2:end)))./(4*pi()^2.*(deltax));
ExactRHS(2,:)= -(tau*deltax.*cos(tau.*xNode(2:end))+sin(tau.*xNode(1:end-1))-sin(tau.*xNode(2:end)))./(4*pi()^2.*(deltax));

%%The LHS requires calculation of \int uh \phi
%which is \int \sum(aN \thetaN) \phi

%Matrix of LHS integrals of form \int \thetaI \phiJ, where \thetaI=\phiI
LHSIntegrals = [1/3 1/6;
                1/6 1/3];
%%Calculate basis weights for approximate solution
%LHSIntegrals * BasisWeights = ExactRHS / deltax
BasisWeights = LHSIntegrals\(ExactRHS./repmat(deltax,size(ExactRHS,1),1));

%%Mapping uh0 values to the appropriate element nodes
uh(1:2:length(xh)) = BasisWeights(1,:);
uh(2:2:length(xh)) = BasisWeights(2,:);

plot(xf,u,'b')
hold on
plot(xh,uh,'r')

BasisWeights_dt = (1./deltax).*(1)