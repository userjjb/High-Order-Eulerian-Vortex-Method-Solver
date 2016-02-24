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
for Ti= 1:6
for Tj= 1:6
fprintf('%i %i | ',Ti,Tj)
Tx=nd2(Ti);
Ty=nd(Tj);

nn=elim(nd',nd',[1 3 2]);
Lag= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,nn(nv,:,:)),bsxfun(@minus,nd(nv),nn(nv,:,:))),3);

K=@(x,y) (y-Ty)./(2*pi*((x-Tx).^2+(y-Ty).^2));

a1= [0,0,Tx,Tx];
a2= [Tx,Tx,1,1];
b1= [0,Ty,0,Ty];
b2= [Ty,1,Ty,1];

for I=1:N
    for J=1:N
        fun=@(x,y) K(x,y).*Lag(x,I).*Lag(y,J); %Lagrange basis for vorticity interp point I,J
        Q=0;
        for Sect=1:4 %Each section with the singularity at the internal corner of all 4
            [X,Y]= meshgrid(linspace(a1(Sect),a2(Sect),it+1),linspace(b1(Sect),b2(Sect),it+1)); %Uniform refinement
            x1=X(1:end-1,1:end-1);  x1=x1(:);
            x2=X(2:end,2:end);      x2=x2(:);
            y1=Y(1:end-1,1:end-1);  y1=y1(:);
            y2=Y(2:end,2:end);      y2=y2(:);

            for it2=1:length(x1)
                if xor(x1>x2,y1>y2) %integral2() will return the wrong sign if the range is reversed
                    C=-1;
                else
                    C=1;
                end
                Q=Q+C*integral2(fun,x1(it2),x2(it2),y1(it2),y2(it2),'RelTol',1e-14,'AbsTol',1e-16);
            end
        end
        Iwmx(J,I,Tj,Ti)=Q;
    end
end
end
fprintf('%i \n',round(toc))
end
save('IwmX6622.mat','IwmX');
warning('on', 'MATLAB:integral2:maxFunEvalsFail')
warning('on', 'MATLAB:integral2:maxFunEvalsPass')
warning('on', 'MATLAB:integral2:minRectSizePass')
warning('on', 'MATLAB:integral2:minRectSizeFail')