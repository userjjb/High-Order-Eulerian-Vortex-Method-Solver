close all
clear all
del = 0.001;
B= 1/2;
A= -1/2;
for N=1:12;
[Qx, Qw]=GLquad(N);
Qx = Qx*(B-A)/2+(B+A)/2;

% xx = A:.0001:B;
% fun = @(x) cos(pi.*x);
% interp = Lagrange(xx,fun,A,B,Qx);
% plot(xx,interp)
I(N)=Qw*cos(pi*Qx);
end
I=I/2;
plot (1:12,log(abs(2/pi-I)))