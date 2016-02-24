clear all
N=5;
f=6/36;%Size of element
f=f/2;

[Qx,Qw] = GLquad(N);
nd=f*(Qx+1);
M=N;
[Qx3,Qw3]= LGLquad(M);
nd3=f*(Qx3+1);
Tx=0.5;%nd3(3);
Ty=0.5;%nd(5);

ndX=nd-Tx;
ndY=nd-Ty;

a=2*(f*0.3*3)^2;
b=2*(f*0.3*3)^2;
dx=f*0.8;
dy=f*0.8;
w=@(x,y) 15*exp(-((x-dx).^2/a+(y-dy).^2/b));%Center is at (dx,dy)

[xx,yy]=meshgrid(nd,nd);
W=ones(numel(nd));%w(xx,yy);

nnX=elim(ndX(1:N)',ndX(1:N)',[1 3 2]);
LagX= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,nnX(nv,:,:)),bsxfun(@minus,ndX(nv),nnX(nv,:,:))),3);
nnY=elim(ndY(1:N)',ndY(1:N)',[1 3 2]);
LagY= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,nnY(nv,:,:)),bsxfun(@minus,ndY(nv),nnY(nv,:,:))),3);
interp_w=@(x,y) reshape(diag(LagX(reshape(x,1,[]),1:N)'*(W'*LagY(reshape(y,1,[]),1:N)))',size(x));

Rsq=@(x,y) (x).^2+(y).^2;
k=@(x,y,del) ((y)./(2*pi*Rsq(x,y))).*(1-(1-(Rsq(x,y)/del^2)).*exp(-Rsq(x,y)/del^2));
del= 0.001*f;
I= integral2(@(x,y) interp_w(x,y).*k(x,y,del),0-Tx,2*f-Tx,0-Ty,2*f-Ty);%,'RelTol',1e-6,'AbsTol',1e-10);

%kk=@(x,y) (y)./(2*pi*((x).^2+(y).^2));
kk=@(x,y) 1./(2*pi*((x./y).*x +y));
kernFun = @(x,y) kk(x,y).*interp_w(x,y);

it3=1;
for it2=2:27
% [a,b]= meshgrid(linspace(0,Tx,it2),linspace(0,Ty,it2)); %Uniform refinement
% [a,b]= meshgrid((Tx)*[1./2.^[0:it2],0],(Ty)*[1./2.^[0:it2],0]); %Tensor semi-refinement
% x1=a(1:end-1,1:end-1);
% x2=a(2:end,2:end);
% y1=b(1:end-1,1:end-1);
% y2=b(2:end,2:end);
%VVV--Quadtree semi-refinement--VVV
seedX= (0+Tx)*[1./4.^[0:it2],0];
seedY= (0+Ty)*[1./4.^[0:it2],0];
x1=[seedX(1:end-1), seedX(2:end-1), seedX(1:end-2)];
x2=[seedX(2:end), zeros(1,it2), seedX(2:end-1)];
y1=[seedY(1:end-1), seedY(1:end-2), seedY(2:end-1)];
y2=[seedY(2:end), seedY(2:end-1), zeros(1,it2)];

Q2=Jint2(x1,x2,y1,y2,xx,yy,f,kernFun,Qw);
sav2(it3)= sum(Q2);
it3=it3+1;
fprintf('%i ',it2)
end
fprintf('\n')




