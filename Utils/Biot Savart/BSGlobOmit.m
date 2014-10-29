function [GlobOmitX, GlobOmit] = BSGlobOmit(fun,A,B,freq,N)
%Global quadrature with singularity omission method for Biot Savart.
%Calculate the velocity at every 1 out of 'freq' quadrature points by 
%integrating the whole domain with a GL quadrature of order N-1 (omitting 
%the eval point in the sum). Where fun is the vorticity density funciton 
%and A and B are the bounds of the self-region.
h = waitbar(0,'Initializing GL points please wait...');
[Qx,Qw]=GLquad(N);
Qx = Qx*(B-A)/2+(B+A)/2;

saved=1;
resolve = 1:round(N/freq):N; %Evaluate at same frequency as Decomp for speed comparison purposes
for i=resolve
    waitbar(saved/numel(resolve),h,'Global Omit method...')
    NOi = [1:i-1,i+1:N];
    GlobOmit(saved) = Qw(NOi)*[ fun(Qx(NOi)) .* (Qx(i)-Qx(NOi))./abs(Qx(i)-Qx(NOi)).^3 ];
    saved=saved+1;
end
GlobOmit = GlobOmit*(B-A)/2;
GlobOmitX = Qx(resolve);
close(h)
end

