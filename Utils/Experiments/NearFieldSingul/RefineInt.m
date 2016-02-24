for it2=2:15
Q=0;
Tx = [-.0625,-.0625,.9375,.9375]; %[-.125,-.125,.875,.875];
Ty = [-.0625,.9375,.9375,-.0625]; %[-.125,.875,.875,-.125];
kk2=@(x,y) (x.*(y-.5))./((x).^2+(y).^2);

for quadr=1:4
% seedX= (Tx(quadr))*[1./1.1.^[0:it2],0];
% seedY= (Ty(quadr))*[1./1.1.^[0:it2],0];
% x1=[seedX(1:end-1), seedX(2:end-1), seedX(1:end-2)];
% x2=[seedX(2:end), zeros(1,it2), seedX(2:end-1)];
% y1=[seedY(1:end-1), seedY(1:end-2), seedY(2:end-1)];
% y2=[seedY(2:end), seedY(2:end-1), zeros(1,it2)];
[a,b]= meshgrid(linspace(0,Tx(quadr),it2),linspace(0,Ty(quadr),it2)); %Uniform refinement
% [a,b]= meshgrid((Tx(quadr))*[1./2.^[0:it2],0],(Ty(quadr))*[1./2.^[0:it2],0]); %Tensor semi-refinement
x1=a(1:end-1,1:end-1); x1=x1(:);
x2=a(2:end,2:end); x2=x2(:);
y1=b(1:end-1,1:end-1); y1=y1(:);
y2=b(2:end,2:end); y2=y2(:);

    for it=1:length(x1)
        if xor(x1>x2,y1>y2)
            C=-1;
        else
            C=1;
        end
        Q=Q+C*integral2(kk2,x1(it),x2(it),y1(it),y2(it),'RelTol',1e-14,'AbsTol',1e-16);
    end
end
end