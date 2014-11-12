clear all
N=1000;
A = magic(N);

%Subtraction speed test
C1 = A-repmat([1:N]',1,N);
%C2 = A-diag(1:N)*ones(size(A)); %VERY SLOW
C3 = bsxfun(@minus,A,[1:N]');
for n=1:N
    C4(n,:) = A-n;
end
%-------------------------------------


%Lagrange weight speed test
[Qx, Qw]=GLquad(N);

NOn = triu(repmat(2:N,N,1))+tril(repmat(1:N-1,N,1),-1);
W1 = prod(bsxfun(@minus,Qx',Qx(NOn')));

for n=1:N
    NOn = [1:n-1 n+1:N];
    W2(n) = prod(Qx(n)-Qx(NOn));
end
min(W2)/eps(min(W2)) %This should be >1e15 or so
%-------------------------------------