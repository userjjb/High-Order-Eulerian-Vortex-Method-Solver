clear all

Np=8;
N=.1*40^2; %Vorticity elements

topm = 4;
Cr = cell(6,topm+1);
Cr{1,1} = 'All operands real';

for m=2:topm+1
if( m == 2 )
    mtimesx('BLAS');
elseif( m == 3 )
    mtimesx('LOOPS');
elseif( m == 4 )
    mtimesx('MATLAB');
elseif( m == 5 )
    mtimesx('SPEED');
elseif( m == 6 )
    mtimesx('LOOPSOMP');
else
    mtimesx('SPEEDOMP');
end
Cr{1,m} = mtimesx;

disp(' ');
disp('--------------------------------------------------------------');
disp('--------------------------------------------------------------');
disp(' ');
disp(['MTIMESX mode: ' mtimesx]);
disp(' ');

disp(' ');
disp('--------------------------------------------------------------');
disp('(real 3x3x1000000) * (real 3x3x1000000) example');
Cr{3,1} = '(3x3xN) *(3x3xN)';
w = rand(Np,Np,N);
Qw = rand(Np,1);
Qwt = rand(1,Np);
r= rand(Np,Np,N,N*Np^2);

% mtimesx
n =4;
tx = zeros(1,n);
for k=1:n
clear Cx
tic
Cx = mtimesx(Qwt,mtimesx(bsxfun(@times,w,r),Qw));
tx(k) = toc;
end
% results
tx = median(tx);

disp(' ');
disp(['MTIMESX Elapsed time ' num2str(tx) ' seconds.']);

end