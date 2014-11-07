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

