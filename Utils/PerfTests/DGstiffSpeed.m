clear all

n=20;
Np=8;
Mp=8;
K= [40 40];

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

disp(['MTIMESX mode: ' mtimesx]);
QwS= rand(Mp,Np,1,Np);
wM= rand(Np,1,Np*K(1)*K(2),1);
vM= rand(1,Mp,Mp*K(1)*K(2),1);

wM0= wM;
wM0(:,1,end/4:end*.75,1)=0;
vM0= vM;
vM0(1,:,end/4:end*.75,1)=0;

w= rand(Np*K(2),Np*K(1));
v= rand(Mp*K(2),Mp*K(1)); 

% mtimesx
tx = zeros(1,n);
for k=1:n
clear Cx
tic
Cx = mtimesx(vM,mtimesx(QwS,wM));
%Cx = mtimesx(vM0,mtimesx(QwS,wM0));
%Cx = mtimesx(vM0(:,:,[1:end/4-1,end*.75+1:end],1),mtimesx(QwS,wM0(:,:,[1:end/4-1,end*.75+1:end],1)));
%Cx= permute(mtimesx(permute(reshape(v',Mp,[]),[3 1 2]) , mtimesx(QwS,permute(reshape(w',Np,[]),[1 3 2])) ),[3 4 1 2]);
tx(k) = toc;
end
% results
tx = median(tx);

disp(['MTIMESX Elapsed time ' num2str(tx) ' seconds.']);

end
L=rand(1,Np);
tic
temp=mtimesx(L,wM);
toc