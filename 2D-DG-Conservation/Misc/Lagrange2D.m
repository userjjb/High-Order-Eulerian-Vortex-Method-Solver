function Interp = Lagrange2D(Reval,fun,A,B,C,D,Qx,Wx,Qy,Wy)
%Evaluates a 2D Lagrange interpolating polynomial for the function 'fun'
%Reval is a E by 2 matrix with E desired evaluation points ( Reval(e,1),Reval(e,2) ) = e_q(x,y)
%Qx is the x-coordinate interpolation points in the tensor grid, Wx being
%the associated barycentric weights (Qy, Wy is the complimentary y-coords)
%With B>A and D>C being the domain edge bounds, x and y coords respectively
%Uses a barycentric evaluation scheme to reuce to O(n) complexity and help
%with high polynomial order float problems (precision)
%Note Qx etc should be ortho poly roots to prevent Runge phenomena (i.e.
%reducing the Lebesgue constant growth)
N= numel(Qx); %Order+1 of interpolations
M= numel(Qy); %^^
nR= size(Reval,1); %Number of eval points

Qx= Qx*(B-A)/2 + (B+A)/2; %Map [-1 1] to [A B]
Qy= Qy*(D-C)/2 + (D+C)/2; %Map [-1 1] to [C D]

zij= fun(reshape(repmat(Qx',M,1),N*M,1),repmat(Qy,N,1)); %Evaluate function to interpolate for each basis 'ij'

%We essentially would like to calculate lX* sum: Wi*zi/deni as a 1D element-wise example
denx= bsxfun(@minus,Reval(:,1),Qx'); %Basis specific denominator that "fixes" lX by dividing it to exclude the unwanted (x-x_i) for basis 'ij'
deny= bsxfun(@minus,Reval(:,2),Qy');
lX= prod(denx,2); %The numerator in the sum
lY= prod(deny,2);
%den=|-------------------denx-------|.*|----deny------|
den= reshape(repmat(denx,M,1),nR,M*N).*repmat(deny,1,N); %Forms total den, we saved denx and deny only so we could calc lX and lY

Wx=reshape(repmat(Wx,M,1),1,N*M); %Form the vectors to calc the tensor product of Wi*Wj
Wy=repmat(Wy,1,N);
%Calculate all eval points elementwise with bsxfun. Not only is this
%compact code and should save on memory (by avoiding needlessly
%constructing the repeated matrices for Wx*Wy), but it is much faster as
%well. ~200 times faster for N=M=50 and size(Reval,1)=100;
Interp= lX.*(bsxfun(@times,Wx,bsxfun(@rdivide,Wy,den))*zij).*lY;
end