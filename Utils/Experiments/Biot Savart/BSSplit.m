function [Split] = BSSplit(fun,A,B,SP,N)
%Split Quadrature method for Biot Savart, for each point SP evaluate the
%velocity by integrating the whole domain in two pieces, the regions left
%and right of the point with a GL quadrature of order N-1. Where fun is the
%vorticity density funciton and A and B are the bounds of the self-region.
h = waitbar(0,'Initializing...');
[Qx,Qw]=GLquad(N);

saved=1;
for i=1:numel(SP) %SP is a set of "split points" to have the velocity evaled at, where the quadrature is divided into left and right portions
    QxL = Qx*(SP(i)-A)/2+(SP(i)+A)/2; %Left set of quadrature points
    QxR = Qx*(B-SP(i))/2+(B+SP(i))/2; %Right set of quadrature points
    waitbar(saved/numel(SP),h,'Split method...')
    Split(1,saved) = (SP(i)-A)/2 * Qw*[ fun(QxL) .* (QxL-SP(i))./abs(SP(i)-QxL).^3 ]; %Eval left integral
    Split(2,saved) = (B-SP(i))/2 * Qw*[ fun(QxR) .* (QxR-SP(i))./abs(SP(i)-QxR).^3 ]; %Eval right integral
    saved=saved+1;
end
close(h)
end

