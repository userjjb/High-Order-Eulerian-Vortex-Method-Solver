function [DecompX, Decomp] = BSDecomp(fun,A,B,K,N)
%Domain decomposition method
[Qx,Qw]=GLquad(N);

Nodes = A:(B-A)/K:B;
El = Nodes(1:end-1);
Er = Nodes(2:end);
deltax = Er-El;
Centroids=(El+Er)/2;

map = Qx*deltax/2+repmat(Centroids,N,1);

saved=1;
resolve = [1:round(max(1,K/800)):K/10,(K/10)+1:round(max(1,K/200)):K-(K/10)-1,K-(K/10):round(max(1,K/800)):K];
%h = waitbar(0,'Initializing...');
for i=resolve
    %waitbar(saved/numel(resolve),h,'Decomposition method...')
    NOi = [1:i-1,i+1:K];
    Decomp(saved) = (Qw*[ fun(map(:,NOi)) .* (Centroids(i)-map(:,NOi))./abs(Centroids(i)-map(:,NOi)).^3 ])*deltax(NOi)'/2;
    saved=saved+1;
end
DecompX = Centroids(resolve);
%close(h)
end