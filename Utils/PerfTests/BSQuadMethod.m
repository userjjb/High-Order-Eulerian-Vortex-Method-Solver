clear all

elems=16^2;
N=5;
[Qx, Qw] = GLquad(N);
w=rand(N);
k=rand(N,N,elems*N*N);

pre=(Qw'*Qw).*w;
preR=pre(:)';
kR=reshape(k,N^2,1,elems*N*N);
for rep=1:2560
    v=mtimesx(Qw,mtimesx(bsxfun(@times,w,k),Qw'));
    v2=sum(sum(bsxfun(@times,pre,k)),2);
    v6=mtimesx(preR,kR);
end