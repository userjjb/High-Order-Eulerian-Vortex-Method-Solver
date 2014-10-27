close all
clear all
clc

%Local vorticity functions
fun = @(x) 1-x.^2/2+x.^4/24-x.^6/720;%+x.^8/40320;
%fun = @(x) 1-(pi^2*x.^2)/2+(pi^4*x.^4)/24-(pi^6*x.^6)/720;  %Taylor trig
%fun = @(x) max(0, exp(-x.^2/.1));                          %Gaussian curve
%fun = @(x) max(0,cos(pi * x)+1 - 2*heaviside(abs(x)-1));   %Cleaved trig
A=-1; %Left boundary
B=1; %Right Boundary

xx=A:(B-A)/200:B;
plot(xx,10*fun(xx),'k') %Exaggerate vorticity function for display purposes

%Domain decomposition method
N=3;
[Qx,Qw]=GLquad(N);

K=100;
Nodes = A:(B-A)/K:B;
El = Nodes(1:end-1);
Er = Nodes(2:end);
deltax = Er-El;
Centroids=(El+Er)/2;

map = Qx*deltax/2+repmat(Centroids,N,1);

saved=1;
resolve = [1:round(max(1,K/800)):K/10,(K/10)+1:round(max(1,K/200)):K-(K/10)-1,K-(K/10):round(max(1,K/800)):K];
h = waitbar(0,'Initializing Decomposition method...');
for i=resolve
    waitbar(saved/numel(resolve),h,'Decomposition method...')
    NOi = [1:i-1,i+1:K];
    Decomp(saved) = (Qw*[ fun(map(:,NOi)) .* (Centroids(i)-map(:,NOi))./abs(Centroids(i)-map(:,NOi)).^3 ])*deltax(NOi)'/2;
    saved=saved+1;
end
hold on
plot(Centroids(resolve),Decomp,'b')
close(h)

%Global quadrature with singularity omission
h = waitbar(0,'Initializing GL points please wait...');
N=500;
[Qx,Qw]=GLquad(N);
Qx = Qx*(B-A)/2+(B+A)/2;
saved=1;
resolve = 1:round(N/numel(Decomp)):N; %Evaluate at same frequency as Decomp for speed comparison purposes
for i=resolve
    waitbar(saved/numel(resolve),h,'Global Omit method...')
    NOi = [1:i-1,i+1:N];
    GlobOmit(saved) = Qw(NOi)*[ fun(Qx(NOi)) .* (Qx(i)-Qx(NOi))./abs(Qx(i)-Qx(NOi)).^3 ];
    saved=saved+1;
end
GlobOmit = GlobOmit*(B-A)/2;
plot(Qx(resolve),GlobOmit,'r')

% %Split Quadrature
% SP = Qx(resolve);
% N=3;
% [Qx,Qw]=GLquad(N);
% saved=1;
% for i=1:numel(SP)
%     QxL = Qx*(SP(i)-A)/2+(SP(i)+A)/2;
%     QxR = Qx*(B-SP(i))/2+(B+SP(i))/2;
%     waitbar(saved/numel(resolve),h,'Split method...')
%     NOi = [1:i-1,i+1:N];
%     Split(1,saved) = (SP(i)-A)/2 * Qw*[ fun(QxL) .* (QxL-SP(i))./abs(SP(i)-QxL).^3 ];
%     Split(2,saved) = (B-SP(i))/2 * Qw*[ fun(QxR) .* (QxR-SP(i))./abs(SP(i)-QxR).^3 ];
%     saved=saved+1;
% end
% plot(SP,Split(1,:)+Split(2,:),'g')

% %Near/Far Split Quadrature
% N=3;
% [Qx,Qw]=GLquad(N);
% saved=1;
% NearFrac = 0.5; %Fraction of split considered near field
% epsil = 0.001;
% for i=1:numel(SP)
%     Lfn=SP(i)-(SP(i)-A)*NearFrac; %Dividing point between near and far fields for left split
%     Rfn=SP(i)+(B-SP(i))*NearFrac; %Dividing point between near and far fields for right split
%     QxLf = Qx*(Lfn-A)/2+(Lfn+A)/2;          %Left far field
%     QxLn = Qx*([SP(i)-epsil]-Lfn)/2+([SP(i)-epsil]+Lfn)/2;  %Left near field
%     QxRn = Qx*(Rfn-[SP(i)+epsil])/2+([SP(i)+epsil]+Rfn)/2;  %Right near field
%     QxRf = Qx*(B-Rfn)/2+(B+Rfn)/2;          %Right far field
%     waitbar(saved/numel(resolve),h,'Near/Far Split method...')
%     NOi = [1:i-1,i+1:N];
%     NFSplit(1,saved) = (Lfn-A)/2 * Qw*[ fun(QxLf) .* (QxLf-SP(i))./abs(SP(i)-QxLf).^3 ];      %Left far field 
%     NFSplit(2,saved) = ([SP(i)-epsil]-Lfn)/2 * Qw*[ fun(QxLn) .* (QxLn-SP(i))./abs(SP(i)-QxLn).^3 ];  %Left near field
%     NFSplit(3,saved) = (Rfn-[SP(i)+epsil])/2 * Qw*[ fun(QxRn) .* (QxRn-SP(i))./abs(SP(i)-QxRn).^3 ];  %Right near field
%     NFSplit(4,saved) = (B-Rfn)/2 * Qw*[ fun(QxRf) .* (QxRf-SP(i))./abs(SP(i)-QxRf).^3 ];      %Right far field
%     saved=saved+1;
% end
% plot(SP,NFSplit(1,:)+NFSplit(2,:)+NFSplit(3,:)+NFSplit(4,:),'c')

close(h)