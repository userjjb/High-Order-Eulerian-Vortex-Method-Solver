function [DecompImp] = BSDecompImp(fun,A,B,VX,deltax,map,Qw)
%Domain decomposition method

saved=1;
for i=1:numel(VX)
    DecompImp(saved) = (Qw*[ fun(map(:,:,i)) .* (VX(i)-map(:,:,i))./abs(VX(i)-map(:,:,i)).^3 ])*deltax(i,:)'/2;
    saved=saved+1;
end

end