function [M,M_] = LagrangeMass(m,Type)
for a=1:m
    for b=1:m
      M(a,b)=InnerProd(m*Type,a,b);
    end
end
M_=inv(M);
end

function q = InnerProd(m,a,b)%,c)
if m>0
    [n,w] = gauss(m);
elseif m<0
    n=[-1,-sqrt(1/5),sqrt(1/5),1];
    m=abs(m);
end


aa=[1:a-1,a+1:m];
bb=[1:b-1,b+1:m];

fun = @(x) ( (x-n(aa(1))).*(x-n(aa(2))).*(x-n(aa(3))) )...
    ./( (n(a)-n(aa(1))).*(n(a)-n(aa(2))).*(n(a)-n(aa(3))) )...
    .*( (x-n(bb(1))).*(x-n(bb(2))).*(x-n(bb(3))) )...
    ./( (n(b)-n(bb(1))).*(n(b)-n(bb(2))).*(n(b)-n(bb(3))) );%...
    %.*( 1./(x-n(c)) );
q= integral(fun,-1,1);
end