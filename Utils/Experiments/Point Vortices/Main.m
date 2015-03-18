clc
close all
clear all

dx=.08;
A=dx^2*10;
[VortX,VortY]=meshgrid(-1:dx:1,-1:dx:1);
a=.25;
b=.25;
VortF=A*exp(-(VortX.^2/a+VortY.^2/b));
mask=VortF>0.02*A;

n=numel(VortF(mask));

VortX=single(reshape(VortX(mask),n,1));
VortY=single(reshape(VortY(mask),n,1));
VortF=single(reshape(VortF(mask),n,1));

%Clone for a vortex pair
VortX=[VortX-4; VortX+4];
VortY=repmat(VortY,2,1);
VortF=[VortF; -VortF];
n=numel(VortF);

figure;
pause(0.00001);
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);
pause(1);

q=1;
dt=0.01;
for t=0:dt:10
    plot(VortX,VortY,'.')
    axis([-6 6 -3 3])
    axis equal
    text(2.02,2,num2str(t));
    pause(0.001)
    
    My=bsxfun(@minus,VortY,repmat(VortY',n,1));
    Mx=bsxfun(@minus,VortX,repmat(VortX',n,1));
    r=(My.^2 + Mx.^2+(dx*q)^2).^(3/2);

    VortX=VortX+(-(My./r)*VortF)*dt;
    VortY=VortY+((Mx./r)*VortF)*dt;
end