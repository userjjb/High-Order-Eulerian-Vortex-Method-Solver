%Evaluate a Nth order Lagrange interpolation for each basis 'n' at points
%'x' returns a n by x matrix: Lx_n(x)
function Lx = Lagrange(N,x,n)
    [Qx Qw] = GLquad(N); %Calc interpolation points
    
    NOn=[1:n-1,n+1:N];
    LagBas_fun = @(x) prod(bsxfun(@minus,x,Qx(NOn))) ./ prod( Qx(n)-Qx(NOn) );
end

%-----------------------------------------
    
%
function [Ld_n] = LagBasisD_n(N,n)

end