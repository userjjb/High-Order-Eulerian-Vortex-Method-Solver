function [ output_args ] = Untitled( input_args )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for i=1:numel(VX)
    frac = min(K-1,max(1,round((VX(i)-del-A)*K/(B-A-2*del)))); %frac out of K elements on left side
    vertices=[A:(VX(i)-del-A)/frac:VX(i)-del , VX(i)+del:(B-VX(i)-del)/(K-frac):B]; %Vertices of all elements, spread evenly for each left and right sides
    centroids = 0.5*( vertices(1:end-1)+vertices(2:end) ); %Centroids of all elements
    KK = find(abs(centroids-VX(i))<eps); %Element that eval point is located in
    noKK = [1:KK-1,KK+1:K+1]; %Remove self element from list
    deltax(i,:) = [ diff(vertices(1:KK)) , diff(vertices(KK+1:end)) ]; %Element widths for all but self
    map(:,:,i) = Qx*deltax(i,:)/2+repmat(centroids(noKK),N,1); %Map GL points onto all elems but self
end

end

