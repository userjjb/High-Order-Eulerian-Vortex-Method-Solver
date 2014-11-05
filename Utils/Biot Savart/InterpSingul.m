close all
clear all
N=450;
[Qx, Qw]=GLquad(N);

del = 0.0001;
xx = del:.001:1;
fun = @(x) 1./(x.^2);
interp = Lagrange(xx,fun,del,1,Qx);
plot(xx,interp)
Qw*(1./(Qx*(1-del)/2+(1+del)/2).^2)

