function [w]=InitialConditions(w,TestCases,wxm,wym)

%2D mollifier with center dx,dy and range over the ellipse with axes a,b
%Not currently used, but may prove useful for IC funs that have 
%discontinuous derivatives at truncation edges for smoothing to avoid the
%interpolation from going haywire
moll=@(x,y,dx,dy,a,b) heaviside(1-(((x-dx)/a).^2+((y-dy)/b).^2)).*exp(1+(-1./(1-(((x-dx)/a).^2+((y-dy)/b).^2).^2)));
%Intial conditions specification-------------------------------------------
ICfuns={}; %Create a cell list of functions that define the ICs
    %Gaussian 1
Ga1=    0.02;
Gb1=    0.04;
Gdx1=   0.12; 
Gdy1=   0.12;
GA1=    .5;
ICfuns{end+1}=@(x,y) GA1*exp(-(((x)-Gdx1).^2/Ga1+(y-Gdy1).^2/Gb1));%Center is at (dx,dy)
    %Gaussian 2
Ga2=    0.012;
Gb2=    0.025;
Gdx2=   -0.12; 
Gdy2=   -0.12;
GA2=    .5;
ICfuns{end+1}=@(x,y) GA2*exp(-(((x)-Gdx2).^2/Ga2+(y-Gdy2).^2/Gb2));
    %Perlmann (3)
Pa=1;
Pb=1;
PR=1;
ICfuns{end+1}=@(x,y) (1/PR^14)*(PR^2-min((x/Pa).^2+(y/Pb).^2,PR^2)).^7;
vAfx=@(x,y) y.*( (-1./(16*(x.^2+y.^2))) .* (1-(1-min(x.^2+y.^2,1)).^8) );
vAfy=@(x,y) -x.*( (-1./(16*(x.^2+y.^2))) .* (1-(1-min(x.^2+y.^2,1)).^8) );
    %Fun (4)
ICfuns{end+1}=@(x,y) moll(x,y,0,0,0.31,0.31);
    %Strain (5-8) patch
S=  [-0.4515, 0.4968, -0.9643, 0.3418];
dx= [-0.6988, 1.4363, -0.1722, -1.5009];
dy= [-1.7756, -1.4566, 0.4175, -0.0937];
p=  [0.6768, 0.3294, 0.5807, 0.2504];
for m=1:4
    ICfuns{end+1}=@(x,y) S(m)*exp(-((x-dx(m)).^2+(y-dy(m)).^2)./p(m)^2);
end
    %Koum \omega^{II} (9)
Ka=1;
Kb=2;
SR=0.8;
ICfuns{end+1}=@(x,y) 1*(1-min( (((x/Ka).^2+(y/Kb).^2).^2)/SR^4 ,1));

%Iterate over each of the IC funs
if TestCases
    for IC=TestCases
        w=w+ICfuns{IC}(wxm,wym);
    end
else
    load('6_36G12_del.mat')
    w=reshape(wxt(:,:,:,279),7*36,7*36)';
end