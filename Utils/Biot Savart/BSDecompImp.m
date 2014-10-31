function [DecompImp] = BSDecompImp(fun,A,B,K,N,VX,del,deltax,map)
%Domain decomposition method
assert(all(VX-del>A) & all(VX+del<B),'Singularity exclusion zone exceeds domain bounds')
%h = waitbar(0,'Initializing...');
[Qx,Qw]=GLquad(N);

% for i=1:numel(VX)
%     frac = min(K-1,max(1,round((VX(i)-del-A)*K/(B-A-2*del)))); %frac out of K elements on left side
%     vertices=[A:(VX(i)-del-A)/frac:VX(i)-del , VX(i)+del:(B-VX(i)-del)/(K-frac):B]; %Vertices of all elements, spread evenly for each left and right sides
%     centroids = 0.5*( vertices(1:end-1)+vertices(2:end) ); %Centroids of all elements
%     KK = find(abs(centroids-VX(i))<eps); %Element that eval point is located in
%     noKK = [1:KK-1,KK+1:K+1]; %Remove self element from list
%     deltax(i,:) = [ diff(vertices(1:KK)) , diff(vertices(KK+1:end)) ]; %Element widths for all but self
%     map(:,:,i) = Qx*deltax(i,:)/2+repmat(centroids(noKK),N,1); %Map GL points onto all elems but self
%     %clc
% end

saved=1;
for i=1:numel(VX)
    %waitbar(saved/numel(VX),h,'Decomposition method...')
    DecompImp(saved) = (Qw*[ fun(map(:,:,i)) .* (VX(i)-map(:,:,i))./abs(VX(i)-map(:,:,i)).^3 ])*deltax(i,:)'/2;
    saved=saved+1;
end
%close(h)
end