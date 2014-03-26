clear all

%%Let u(x,t) be the exact solution on the domain 0<=x<=1
%%Let u(.,0) = u0 = sin(2pi * x)

%%Let vh be the finite vector space of linear polynomials
%%Let uh be the approximate numerical solution consisting of a linear 
%combination of basis functions (\thetaN for the Nth basis) in vh with
%scalar coefficients aN (i.e. \sum aN \thetaN)
%%Let \phiN be the Nth test function in the same vector space as the basis
%functions (vh)

%%Discretize the domain into K elements with K+1 nodes
tau=2*pi();

K=10;
x=0:1/K:1;
deltax = diff(x);
%All the nodes, note their are repeats since nodes have a coincident
%brother from an adjacent element, except at the boundaries
xh = sort([x,x(2:end-1)]);

xu = 0:1/100:1;
u= sin(tau.*xf);

%%According to Cockburn,Shu 2001 Eqn 2.2 let uh(.,0) be computed by
%\int uh \phi = \int u0 \phi for each element (xj-1/2 < x < xj=1/2)
%%We can explicitly define a formula for cN = \int u0 \phiN the value of
%the RHS
%%NOTE: There is an extra factor of 1/deltax here that was moved from the LHS
c1= (tau*deltax.*cos(tau.*x(1:end-1))+sin(tau.*x(1:end-1))-sin(tau.*x(2:end)))./(4*pi()^2.*(deltax).^2);
c2= -(tau*deltax.*cos(tau.*x(2:end))+sin(tau.*x(1:end-1))-sin(tau.*x(2:end)))./(4*pi()^2.*(deltax).^2);

%%The LHS requires \int uh \phi which is \int \sum(aN \thetaN) \phi
%Since both \phi and \theta belong to vh this becomes:
%\sum aN \int \theta_N \theta_M
%where possibly M=N.
%%Here N<=2 so we only need two aN for the two possibilities for the 
%integral.
%%We find that the LHS integral is either h/3 for M=N or h/6 for M!=N, we
%%then solve algebraically for aN
a1 = 4*c1 - 2*c2;
a2 = 4*c2 - 2*c1;

%%Mapping uh0 values to the appropriate element nodes
uh(1:2:length(xh)) = a1;
uh(2:2:length(xh)) = a2;

plot(xu,u)
hold on
plot(xh,uh)