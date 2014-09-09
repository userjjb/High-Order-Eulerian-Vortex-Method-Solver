function [M,M_] = LagrangeMass(m,Type)
switch Type
    case 'GL'
        [n,~] = GL_Quad(m);
    case 'LGL'
        [n,~,~] = LGL_Quad(m);
    otherwise
        return
end

for a=1:m
    for b=1:m
      M(a,b)=InnerProd(n,a,b);
    end
end
M_=inv(M);
end

function q = InnerProd(n,a,b)
m=numel(n);

aa=[1:a-1,a+1:m];
bb=[1:b-1,b+1:m];

fun = @(x) prod(bsxfun(@rdivide, bsxfun(@minus,x,n(aa)) ,n(a)-n(aa)))...
    .*prod(bsxfun(@rdivide, bsxfun(@minus,x,n(bb)) ,n(b)-n(bb)));

q= integral(fun,-1,1);
end