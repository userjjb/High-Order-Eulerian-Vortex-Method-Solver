%Investigates behavior of polynomial interpolations of fluxes where the
%velocity is non-constant across the element, e.g. f(x)=c(x)*u(x)
%
%One could interpolate c and u separately and multiply them together to
%arrive at an interpolation of c(x)*u(x), but if each individual interpolation is
%of order n, then the resulting product is of order 2*n
%
%Alternatively one could interpolate the value of the product f(x) with
%only order n.
%
%Hypothesis: The interpolation of the product will converge faster than the
%product of the interpolations.
%Reasoning: It is a weaker condition to specify the flux interpolation as
%the result of a product of interpolations. It ends up being a "happy
%coincidence" that c(x)*u(x)=c_h(xi)*u_h(xi) rather than an explicitly
%enforced condition
%
%Added benefit: The integral of the flux in a FEM-style stiffness matrix 
%can be approximated with a Gaussian quadrature of the flux interpolation
%times the spatial derivative of the shape function. If the interpolation
%of the product is used the interpolation points can be collocated with the
%quadrature points, simplifying the sums.
%
%For our test velocity and conserved quantity functions we use a sine
%function and Gaussian bell curve to try and provide a challenging
%interpolation. In reality the variation of velocity inside of an element
%should be much better well behaved and a lower order interpolation would
%be sufficient to reach machine precision compared to the test here
clc
close all
%Discretization and exact functions
x = -1:1/500:1;
f = sin(pi*x);
c2 = 20/100;
b=-.2;
g = exp(-(x-b).^2/(2*c2));
f_g=f.*g;

%------------------
%Order m-1 interpolating polynomial, using Gauss-Legendre nodes
m=11;
[nd,w] = gauss(m);

%Interpolation of each example function using a Lagrange basis
f_h=zeros(1,1001);
for a=1:m
    aa=[1:a-1,a+1:m];
    multer=1;
    for i=1:m-1
        %Calculate a Lagrange basis by the product of terms (x-xi)/(xj-xi)
        multer= multer.* ( ( x-nd(aa(i)) )/( nd(a)-nd(aa(i)) ) );
    end
    %Assemble intepolation values based on function values at nodes
    f_h = f_h + sin(pi*nd(a)).*multer;
end

g_h=zeros(1,1001);
for a=1:m
    aa=[1:a-1,a+1:m];
    multer=1;
    for i=1:m-1
        multer= multer.* ( ( x-nd(aa(i)) )/( nd(a)-nd(aa(i)) ) );
    end
    g_h = g_h + exp(-(nd(a)-b).^2/(2*c2)).*multer;
end

%Order 2(m-1) order interpolation resulting from the product of
%interpolants
fh_gh = f_h.*g_h;

%Order m-1 interpolation of the product
f_g_h=zeros(1,1001);
for a=1:m
    aa=[1:a-1,a+1:m];
    multer=1;
    for i=1:m-1
        multer= multer.* ( ( x-nd(aa(i)) )/( nd(a)-nd(aa(i)) ) );
    end
    f_g_h = f_g_h + sin(pi*nd(a)).*exp(-(nd(a)-b).^2/(2*c2)).*multer;
end

plot(x,f_g,'b')
hold on
plot(x,f_g_h,'r')
plot(x,fh_gh,'g')

%Report norms for comparison
norm(f_g-fh_gh)
norm(f_g-f_g_h)
