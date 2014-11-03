function [NFSplit] = BSNFSplit(fun,A,B,SP,N,NearFrac,epsil)
%Near/Far Split Quadrature
[Qx,Qw]=GLquad(N);

saved=1;
for i=1:numel(SP)
    Lfn=SP(i)-(SP(i)-A)*NearFrac; %Dividing point between near and far fields for left split
    Rfn=SP(i)+(B-SP(i))*NearFrac; %Dividing point between near and far fields for right split
    QxLf = Qx*(Lfn-A)/2+(Lfn+A)/2;          %Left far field
    QxLn = Qx*([SP(i)-epsil]-Lfn)/2+([SP(i)-epsil]+Lfn)/2;  %Left near field
    QxRn = Qx*(Rfn-[SP(i)+epsil])/2+([SP(i)+epsil]+Rfn)/2;  %Right near field
    QxRf = Qx*(B-Rfn)/2+(B+Rfn)/2;          %Right far field
    NOi = [1:i-1,i+1:N];
    NFSplit(1,saved) = (Lfn-A)/2 * Qw*[ fun(QxLf) .* (QxLf-SP(i))./abs(SP(i)-QxLf).^3 ];      %Left far field 
    NFSplit(2,saved) = ([SP(i)-epsil]-Lfn)/2 * Qw*[ fun(QxLn) .* (QxLn-SP(i))./abs(SP(i)-QxLn).^3 ];  %Left near field
    NFSplit(3,saved) = (Rfn-[SP(i)+epsil])/2 * Qw*[ fun(QxRn) .* (QxRn-SP(i))./abs(SP(i)-QxRn).^3 ];  %Right near field
    NFSplit(4,saved) = (B-Rfn)/2 * Qw*[ fun(QxRf) .* (QxRf-SP(i))./abs(SP(i)-QxRf).^3 ];      %Right far field
    saved=saved+1;
end
end

