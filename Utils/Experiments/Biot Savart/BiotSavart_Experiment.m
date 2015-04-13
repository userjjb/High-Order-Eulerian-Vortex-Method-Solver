close all
clear all
clc

%Local vorticity functions
%fun = @(x) 1-x.^2/2+x.^4/24-x.^6/720+x.^8/40320-x.^10/3628800; %Taylor trig
%fun = @(x) 1-(pi^2*x.^2)/2+(pi^4*x.^4)/24-(pi^6*x.^6)/720;     %Taylor trig
%fun = @(x) max(0, exp(-x.^2/.04));                              %Gaussian curve
fun = @(x) max(0,cos(pi*x)+1- heaviside(abs(x)-1.01));       %Cleaved trig
%fun = @(x) heaviside(x+0.5).*exp(-1./(1-(2*x).^2)).*heaviside(0.5-x); %Mollifier
A=-1; %Left boundary
B=1; %Right Boundary

xx=A:(B-A)/200:B;
hold on
plot(xx,10*fun(xx),'k') %Exaggerate vorticity function for display purposes

%Domain decomposition method
[DecompX, DecompV] = BSDecomp(fun,A,B,200,3);
plot(DecompX,DecompV,'b')

%Improved decomp method with hp-adaptivity
del2 = (B-A)/1000;
truncated = DecompX(DecompX-del2>A&DecompX+del2<B);
[deltax, map, Qw] = AdapElem(A,B,100,10,del2,truncated);
[DecompImp] = BSDecompImp(fun,A,B,truncated,deltax,map,Qw);
plot(truncated,DecompImp,'m')
text(truncated(find(max(DecompImp)==DecompImp)),max(DecompImp),num2str(max(DecompImp)));

% %Correct for external vorticity for each of the resolved points
% cut=1.048; %How far into external domain corrected quadratures extend
% N=100;
% [Qx,Qw]=GLquad(N);
% QxR = Qx*(pi/cut-B)/2+(pi/cut+B)/2;
% QxL = Qx*(A--pi/cut)/2+(A+-pi/cut)/2;
% saved=1;
% for i=Centroids(resolve)
%     ExtR(saved)=Qw*[ fun(QxR) .* (i-QxR)./abs(i-QxR).^3 ];
%     ExtL(saved)=Qw*[ fun(QxL) .* (i-QxL)./abs(i-QxL).^3 ];
%     saved=saved+1;
% end
% plot(Centroids(resolve),Decomp+ExtR+ExtL,':b')

N=200;
del = (B-A)/200;
[Qx,Qw]=GLquad(N);
%Global quadrature with singularity omission
[GlobOmitX, GlobOmit] = BSGlobOmit(fun,A,B,numel(DecompV),N,Qx,Qw);
plot(GlobOmitX,GlobOmit,'r')

[RoseMooreX, RoseMoore] = BSRoseMoore(fun,A,B,numel(DecompV),N,Qx,Qw,del);
plot(RoseMooreX,RoseMoore,'g')

[WinLeonX, WinLeon] = BSWinLeon(fun,A,B,numel(DecompV),N,Qx,Qw,del);
plot(WinLeonX,WinLeon,'c')

%axis([0,1,0,50])

%Split Quadrature
% [Split] = BSSplit(fun,A,B,GlobOmitX,6);
% plot(GlobOmitX,Split(1,:)+Split(2,:),'--g')

% %Near/Far Split Quadrature method
% [NFSplit] = BSNFSplit(fun,A,B,GlobOmitX,100,0.01,0);
% plot(GlobOmitX,NFSplit(1,:)+NFSplit(2,:)+NFSplit(3,:)+NFSplit(4,:),':b')


