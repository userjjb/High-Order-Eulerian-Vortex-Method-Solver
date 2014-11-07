function Interp = Lagrange2D(Reval,fun,A,B,C,D,Qx,Wx,Qy,Wy)
%Reval is a Q by 2 matrix with Q evaluation points
Qx = Qx*(B-A)/2 + (B+A)/2; %Map [-1 1] to [A B]
Qy = Qy*(D-C)/2 + (D+C)/2; %Map [-1 1] to [A B]
nQx= numel(Qx);
nQy= numel(Qy);
zij = fun(reshape(repmat(Qx',nQy,1),nQx*nQy,1),repmat(Qy,nQx,1)); %Evaluate function to interpolate for each basis 'ij'
nR= size(Reval,1);

denx = bsxfun(@minus,Reval(:,1),Qx'); %Denominator
deny = bsxfun(@minus,Reval(:,2),Qy');
lX = prod(denx,2);
lY = prod(deny,2);
%den=|-------------------denx-------------|.*|----deny---------|
den= reshape(repmat(denx,nQy,1),nR,nQy*nQx).*repmat(deny,1,nQx);

nWx= numel(Wx);
nWy= numel(Wy);
Wx=reshape(repmat(Wx,nWy,1),1,nWx*nWy);
Wy=repmat(Wy,1,nWx);
Interp = lX.*(bsxfun(@times,Wx,bsxfun(@rdivide,Wy,den))*zij).*lY;
end