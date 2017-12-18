N=6;
[Qx,Qw] = GLquad(N);
[Qx2,Qw2]= LGLquad(N);
nd=(Qx+1)/2;
nd2=(Qx2+1)/2;

for Ex=1:3
for Ey=1:3

for Ti= 1:6
for Tj= 1:6

Tx=nd2(Ti);
Ty=nd(Tj);

K=@(x,y) ((y+(Ey-2))-Ty)./(2*pi*(((x+(Ex-2))-Tx).^2+((y+(Ey-2))-Ty).^2));

for I=1:N
    for J=1:N
        KK(J,I,Tj,Ti,Ex,Ey)=K(nd(J),nd(I));
    end
end

end
end
end
end

KM = -4*bsxfun(@rdivide,Iwm(:,:,3,3,:,:,2),Qw'*Qw);
KKK = KK(:,:,3,3,:,:);

D=KM-KKK;

for Ey=1:3
    for Ex=1:3
        KKK2((Ey-1)*6+1:(Ey-1)*6+6,(Ex-1)*6+1:(Ex-1)*6+6)=KKK(:,:,:,:,Ey,Ex);
        KM2((Ey-1)*6+1:(Ey-1)*6+6,(Ex-1)*6+1:(Ex-1)*6+6)=KM(:,:,:,:,4-Ex,4-Ey); %Reverse order compared to KKK and sqapped Ex,Ey
    end
end

nodes=bsxfun(@plus,Qx/2,[-1,0,1]);
[xx,yy] = meshgrid(nodes(:),nodes(:));
figure(1)
surf(xx,yy,-log10(abs(KKK2-KM2)),'facecolor', 'interp')
view(2);
axis equal;
cmap=viridis();
colormap(flipud(cmap))
axis([-1.5,1.5,-1.5,1.5])
colorbar()

figure(2)
surf(xx,yy,KKK2,'facecolor', 'interp')
view(2);
axis equal;
cmap=viridis();
colormap(flipud(cmap))
axis([-1.5,1.5,-1.5,1.5])
colorbar()