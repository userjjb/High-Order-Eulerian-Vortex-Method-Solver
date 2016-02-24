close all; clear all; clc
warning('off', 'MATLAB:integral2:maxFunEvalsFail')
warning('off', 'MATLAB:integral2:maxFunEvalsPass')
warning('off', 'MATLAB:integral2:minRectSizePass')
warning('off', 'MATLAB:integral2:minRectSizeFail')
it=1; %May want this to be more in the future (or a loop over a range of 'it', though it looks 
        %like this is ineffective at least without better heuristics to determine the "best" result)

N=6;
[Qx,Qw] = GLquad(N);
[Qx2,Qw2]= LGLquad(N);
nd=(Qx+1)/2;
nd2=(Qx2+1)/2;
tic
for Ex=1:3
for Ey=1:3
    fprintf('%i %i\n',Ex,Ey)
for Ti= 1:6
for Tj= 1:6
fprintf('%i %i | ',Ti,Tj)
Tx=nd(Ti);
Ty=nd2(Tj);

nn=elim(nd',nd',[1 3 2]);
Lag= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,nn(nv,:,:)),bsxfun(@minus,nd(nv),nn(nv,:,:))),3);

K=@(x,y) ((x+(Ex-2))-Tx)./(2*pi*(((x+(Ex-2))-Tx).^2+((y+(Ey-2))-Ty).^2));

for I=1:N
    for J=1:N
        fun=@(x,y) K(x,y).*Lag(x,I).*Lag(y,J); %Lagrange basis for vorticity interp point I,J
        IwmyN(J,I,Tj,Ti,Ex,Ey)=integral2(fun,0,1,0,1,'RelTol',1e-14,'AbsTol',1e-16);
    end
end
end
fprintf('%i \n',round(toc))
end
end
end
save('IwmYN6622.mat','IwmyN');
warning('on', 'MATLAB:integral2:maxFunEvalsFail')
warning('on', 'MATLAB:integral2:maxFunEvalsPass')
warning('on', 'MATLAB:integral2:minRectSizePass')
warning('on', 'MATLAB:integral2:minRectSizeFail')