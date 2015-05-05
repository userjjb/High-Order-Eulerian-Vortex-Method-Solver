[X,Y] = meshgrid (0:.1:2*pi,0:.1:2*pi);

dt = 0.01;
endt = 5;

c=2*pi;
theta = pi/4;
cx = c*cos(theta);
cy = c*sin(theta);

h=figure;set(gcf, 'Color','white')
nFrames = length(endt/dt+1);
vidObj = VideoWriter('2DConvDG_B.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 1/(2*dt);
open(vidObj);

for t=0:dt:endt
    Z = sin(X-cx*t).*sin(Y-cy*t);
    surf(X,Y,Z)
    view([-9 67])
    text(0,8,1.2,'Time:')
    text(0.65,8,1.11,num2str(t));
    writeVideo(vidObj, getframe(h));
    %pause(0.0001)
end

close(gcf)
close(vidObj);