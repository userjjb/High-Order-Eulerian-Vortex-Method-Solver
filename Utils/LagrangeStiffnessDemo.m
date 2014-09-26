function saved = LagrangeStiffness(n)
saved=zeros(n,2);
for q=1:n
    [exact,interp] = Scratch(n,q);
    saved(q,1) = exact;
    saved(q,2) = interp;
end
res = norm(saved(:,1)-saved(:,2))
end

%-----------------------------------------
    
%Calculate stiffness
function [exact,interp,Ld_q] = Scratch(m,q)
[nd,w]=gauss(m);

c2 = 20/100;
b=-.2;

NOq=[1:q-1,q+1:m];
for i=1:m
    NOqi=NOq(not(NOq==i));
    dLagBas_fun = @(x) prod( bsxfun(@minus,x,nd(NOqi)) )./prod( nd(q)-nd(NOq) );
    if not(i==q)
        Ld_q(i) = dLagBas_fun(nd(i));
    elseif i==q
        Ld_q(i) = sum( 1./( nd(i)-nd(NOq) ) );
    end
end

i=1:m;
interp = sum( Ld_q(i)'.*w(i)'.*sin(pi*nd(i)).*exp(-(nd(i)-b).^2/(2*c2)) );

f_g_ld_fun = @(x) sin(pi.*x).*exp(-(x-b).^2./(2*c2)).*prod(bsxfun(@rdivide, bsxfun(@minus,x,nd(NOq)) ,nd(q)-nd(NOq))).*sum(1./bsxfun(@minus,x,nd(NOq)) );
exact= integral(f_g_ld_fun,-1,1);
end