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
clear all
format long
%Exact functions
s = @(x) sin(pi*x);

c2 = 20/100;
b=-.2;
g = @(x) exp(-(x-b).^2/(2*c2));

s_g=@(x) s(x).*g(x);

%Quad stuff
N=11;
[nd,w] = gauss(N);

%Build out the funs for our Lagrange basis and it's derivative
Lag= @(x,n) prod(bsxfun(@rdivide,bsxfun(@minus,x,nd([1:n-1,n+1:N])),bsxfun(@minus,nd(n),nd([1:n-1,n+1:N]))));
dift= @(x,n) sum(1./bsxfun(@minus,x,nd([1:n-1,n+1:N])));
dLag= @(x,n) Lag(x,n).*dift(x,n);

%I haven't come up with a clever way to generate an arbitrary order interp
%fun on the fly, this will have to do...
%many term matrix saved here for convenience, delete terms as needed
%ev=@(x) [Lag(x,1)',Lag(x,2)',Lag(x,3)',Lag(x,4)',Lag(x,5)',Lag(x,6)',Lag(x,7)',Lag(x,8)',Lag(x,9)',Lag(x,10)',Lag(x,11)',Lag(x,12)',Lag(x,13)',Lag(x,14)',Lag(x,15)',Lag(x,16)']';
ev=@(x) [Lag(x,1)',Lag(x,2)',Lag(x,3)',Lag(x,4)',Lag(x,5)',Lag(x,6)',Lag(x,7)',Lag(x,8)',Lag(x,9)',Lag(x,10)',Lag(x,11)']';

interp_s=@(x) s(nd')*ev(x);
interp_g=@(x) g(nd')*ev(x);
interp_sg=@(x) interp_s(x).*interp_g(x);    %Product of interps

interp_s_g=@(x) s_g(nd')*ev(x);             %Interp of product

n=6; %Particular Lagrange basis function n<=N
stiff= @(x) s_g(x).*dLag(x,n);              %Exact-ish term to be used in stiffness integral
stiff_isg= @(x) interp_sg(x).*dLag(x,n);    %Product of interps
stiff_is_g= @(x) interp_s_g(x).*dLag(x,n);  %Interp of product

%Integral of flux
%Let's confirm that the two methods should be equal here(interestingly IMO)
%and not just an artifact of quadrature
fprintf('Exact-ish Integral, followed by deviation of methods\n');
fprintf('%.17f \n',-0.1234568901234567890);%Dummy number for precision header
res(1)=integral(s_g,-1,1,'RelTol',1e-14,'AbsTol',1e-16);        %Exact-ish
fprintf('%.17f Integral of SG\n',res(1));
res(2)=integral(interp_sg,-1,1,'RelTol',1e-14,'AbsTol',1e-16);    %Product of interps
fprintf('%.17f Integral of Prod of In\n',res(1)-res(2));
res(3)=integral(interp_s_g,-1,1,'RelTol',1e-14,'AbsTol',1e-16);   %Interp of product
fprintf('%.17f Integral of In of Prod\n',res(1)-res(3));
%They are equal! And not just because of the quadrature, cool!

%Stiffness entries
%Here be the dragons, ideally our Matlab integrals would be equal, and
%close to the quadrature
fprintf('Exact-ish Stiffness, followed by deviation of methods\n');
res(4)=integral(stiff,-1,1,'RelTol',1e-14,'AbsTol',1e-16);        %Exact-ish
fprintf('%.17f Stiffness\n',res(4));
res(5)=integral(stiff_isg,-1,1,'RelTol',1e-14,'AbsTol',1e-16);    %Product of interps
fprintf('%.17f Stiff: Prod of In\n',res(4)-res(5));
res(6)=integral(stiff_is_g,-1,1,'RelTol',1e-14,'AbsTol',1e-16);   %Interp of product (worse than ^)
fprintf('%.17f Stiff: In of Prod\n',res(4)-res(6));
res(7)=( s(nd').*g(nd').*dLag(nd'-eps(nd'),n) )*w'; %Equal to interp of product
fprintf('%.17f Stiff: In of Prod (Quad)\n',res(4)-res(7));
%Result: The Matlab integrals aren't equal. This is a consequence of the
%fact that although \int L_i = \int L_i L_j
%it is NOT equal to \int L_i L_j L_k, the result of trying to perform
%quadrature on the triple product of three interps

%Examine how well the interp fits the exact
xx=-1:.01:1;
plot(xx,stiff(xx),'k')
hold on
plot(xx,stiff_isg(xx),'g')
plot(xx,stiff_is_g(xx),'r')

%What if we consider the effect the quad of three interps has?
%At least here we hope we can recover the Matlab integral of the product of
%the interps
for i=1:N
    for j=1:N
        for k=1:N
            trip=@(x) Lag(x,i).*Lag(x,j).*Lag(x,k);
            W(i,j,k)=integral(trip,-1,1,'RelTol',1e-14,'AbsTol',1e-16);
        end
    end
end
tot1=0;
for i=1:N
    for j=1:N
        for k=1:N
            tot1=tot1+( s(nd(i))*g(nd(j))*dLag(nd(k)-eps(1),6)*W(i,j,k) );
        end
    end
end
res(8)=tot1; %As expected: equal to product of interps
fprintf('%.17f Stiff:Prod of In (Quad N^3)\n',res(4)-res(8));
%It sucks that we move from needing N weights, to N^3 weights: W(i,j,k)
%We can reduce comp. effort for initialization by about 80% by only 
%computing unique triplets, that is \int(1,2,3)=\int(3,2,1) etc.

%What if we include the Ln(x)' in our quad weights?
for i=1:N
    for j=1:N
            trip=@(x) Lag(x,i).*Lag(x,j).*dLag(x,n);
            W2(i,j)=integral(trip,-1,1,'RelTol',1e-14,'AbsTol',1e-15);
    end
end
tot2a=0;
%This is the slower way to eval with loops
for i=1:N
    for j=1:N
            tot2a=tot2a+( s(nd(i))*g(nd(j))*W2(i,j) );
    end
end
res(9)=tot2a;
fprintf('%.17f Stiff:Prod of In (Quad N^2 for() )\n',res(4)-res(9));
%We can vectorize this to
res(10)=(s(nd')*W2)*g(nd);
fprintf('%.17f Stiff:Prod of In (Quad N^2 vect)\n',res(4)-res(10));
%Equal to tot1, no better or worse
%We reducted it to N^2 weights: W(i,j),but need one for each dL_n, back to
%the N^3 required from above
%We can reduce comp. effort for initialization by about 50% by only 
%computing unique doublets, that is \int(1,2,n)=\int(2,1,n) etc.
%Not as cheap initialization as tot1 method, but we don't need to store 
%dLag(x_i) and our W2(:,:,n) is symmetric. It is also much more convenient
%to calculate the complete sum

%How good can we do if we use our own quad to generate W3~=W2, instead of
%relying on Matlab 'integral()'
[nd2,w2] = gauss(N+1);
for i=1:N
    for j=1:N
            W3(i,j)=(Lag(nd2',i).*Lag(nd2',j).*dLag(nd2',n))*w2';
    end
end
res(11)=(s(nd')*W3)*g(nd);
fprintf('%.17f Stiff:Prod of In MyQuad(Quad N^2 vect)\n',res(4)-res(11));

%Exact polynomial expressions for N=11,n=6 and N=9,n=5 for validating
%against. Rather than functions these should be done as symfuns to be more
%exact. These take advantage of symmetric basis functions about x=0 to cut
%the needed terms in half (otherwise I'd pull my hair out). A=Qx(1)^2,
%B=Qx(2)^2,...
%We can also auto-generate something like these with our binomial expansion
%program in PolyVandLag.m, but it suffers from accumulated precision errors
%Lag11_6= @(x,A,B,C,D,E) (x.^10 - x.^8*(A+B+C+D+E) + x.^6*(A*B + (A+B)*C + (A+B+C)*D + (A+B+C+D)*E) - x.^4*(D*E*(B+C) + B*C*(D+E) + A*(B*(C+D+E)+C*(D+E)+D*E)) + x.^2*(C*D*E*(A+B) + A*(B*(C*(D+E)+D*E))) - A*B*C*D*E)/-(A*B*C*D*E);
%Lag9_5= @(x,A,B,C,D) (x.^8 - x.^6*(A+B+C+D) + x.^4*(A*B+C*(A+B)+D*(A+B+C)) - x.^2*(A*B*C+D*(A*B+C*(A+B))) + A*B*C*D)/-(A*B*C*D)