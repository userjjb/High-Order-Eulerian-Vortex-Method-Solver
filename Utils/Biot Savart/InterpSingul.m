close all
clear all
xx=-1:.01:1;
N=2000;
[Qx, ~]=GLquad(N);
% N=34
% del=1;
% 
% [Qx,Qw] = GLquad(N);
% LEval = Lagrange(N,xx,1:N);
% b = (Qx').^2;
% 
% interp=b*LEval;
% plot(xx,interp)
% Qw*(1./(1+del-Qx).^2)

for n=1:N
    NOn = [1:n-1 n+1:N];
    den = prod(Qx(n)-Qx(NOn));
    iter=1;
    for i=xx
        val(n,iter)=prod(i-Qx(NOn))/den;
        iter=iter+1;
    end
end