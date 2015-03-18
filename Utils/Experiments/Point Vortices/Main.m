clc
close all
clear all

dx=.06;
A=dx^2*10;
[VortX,VortY]=meshgrid(-2:dx:2,-2:dx:2);
a=.5;
b=1;
VortF=A*exp(-(VortX.^2/a+VortY.^2/b));
mask=VortF>0.01*A;

n=numel(VortF(mask));

VortX=single(reshape(VortX(mask),n,1));
VortY=single(reshape(VortY(mask),n,1));
VortF=single(reshape(VortF(mask),n,1));

figure;
pause(0.00001);
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);
pause(1);

q=1;
dt=0.005;
for t=0:dt:1
    plot(VortX,VortY,'.')
    axis([-3 3 -3 3])
    axis square
    text(2.02,2,num2str(t));
    pause(0.001)
    
    My=bsxfun(@minus,VortY,repmat(VortY',n,1));
    Mx=bsxfun(@minus,VortX,repmat(VortX',n,1));
    r=(My.^2 + Mx.^2+(dx*q)^2).^(3/2);

    VortX=VortX+(-(My./r)*VortF)*dt;
    VortY=VortY+((Mx./r)*VortF)*dt;
end