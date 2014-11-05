function [interp, W] = Lagrange(Xeval,fun,A,B,Qx,W)
%Evaluates the Nth order Lagrange interpolating polynomial at all points 
%'Xeval' for the function handle 'fun' in the domain [A,B]
%Xeval must not equal any Qx
%It is a good idea to save W for future use, if W is left out of the input
%arguments it will be calculated.
%Able to operate in O(n) time by using barycentric evaluation
%Qx should be the Legendre root points (to avoid problems of ill-conditioning)
%N (for GLquad(N)) must be <765, otherwise precision errors corrupt the result
Qx = Qx*(B-A)/2 + (B+A)/2; %Map [-1 1] to [A B]
yj = fun(Qx); %Evaluate function to interpolate for each basis 'j'

if nargin==5
    W = LagBaryWeight(Qx); %Calculate barycentric weights, reuse in future
end
den = bsxfun(@minus,Xeval,Qx); %Denominator
lX = prod(den);

interp = lX'.*((repmat(W,numel(Xeval),1)./den')*yj);
end

function W = LagBaryWeight(Qx)
%Generates the barycentric weights for an interpolating Lagrange polynomial
%with interpoaltion points Qx (typically ortho poly roots)
    N=numel(Qx);
    assert(N<1000,'Order of Lagrange polynomial must be <1000 to avoid precision errors')
    for n=1:N
        NOn = [1:n-1 n+1:N];
        W(n) = 1/prod(Qx(n)-Qx(NOn));
    end
end