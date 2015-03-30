clear all

n=4;
a=4; %Np
c=a*100*100; %Np x elems_x x elems_y
nn=a; %W_Np

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
A = rand(a,a,c,nn);
B = rand(a,1,c);
BB = rand(1,a,c);

% mtimesx
tx = zeros(1,n);
for k=1:n
clear Cx
tic
Cx = mtimesx(BB,mtimesx(A,B));
tx(k) = toc;
end
% results
tx = median(tx);

disp(' ');
disp(['MTIMESX Elapsed time ' num2str(tx) ' seconds.']);

end