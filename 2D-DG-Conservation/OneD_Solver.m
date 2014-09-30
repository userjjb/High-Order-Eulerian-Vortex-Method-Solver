N =5; %Polynomial order
Np =N+1; %Number of quadrature and interpolation points

K =11; %Number of elements
Nd = 2*pi*(0:1/(K-1):1); %Element nodes

[Qx,Qw] = GLquad(N); %Gauss Legendre quadrature points and weights

Stiffness = LagrangeStiffness(N);