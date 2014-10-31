function [RoseMooreX, RoseMoore] = BSRoseMoore(fun,A,B,freq,N,del)
%Global quadrature with singularity omission method for Biot Savart.
%Calculate the velocity at every 1 out of 'freq' quadrature points by 
%integrating the whole domain with a GL quadrature of order N-1 (omitting 
%the eval point in the sum). Where fun is the vorticity density funciton 
%and A and B are the bounds of the self-region.
%h = waitbar(0,'Initializing GL points please wait...');
[Qx,Qw]=GLquad(N);
Qx = Qx*(B-A)/2+(B+A)/2;

saved=1;
resolve = 1:max(1,round(N/freq)):N; %Evaluate at same frequency as Decomp for speed comparison purposes
for i=resolve
    %waitbar(saved/numel(resolve),h,'Rosenhead-Moore method...')
    %NOi = [1:i-1,i+1:N];
    RoseMoore(saved) = Qw*[ fun(Qx) .* (Qx(i)-Qx)./(abs(Qx(i)-Qx).^2+del^2).^(3/2) ];
    saved=saved+1;
end
RoseMoore = RoseMoore*(B-A)/2;
RoseMooreX = Qx(resolve);
%close(h)
end

