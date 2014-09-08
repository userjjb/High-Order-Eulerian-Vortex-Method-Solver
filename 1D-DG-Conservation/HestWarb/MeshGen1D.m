function [Nv, VX, K, EToV] = JMeshGen1D(xmin,xmax,K)

% function [Nv, VX, K, EToV] = MeshGen1D(xmin,xmax,K)
% Generate simple equidistant grid with K elements

%Number of vertices
Nv = K+1; 

% Generate node coordinates
VX = (xmax-xmin)*(0:Nv-1)/(Nv-1) + xmin;

% read element to node connectivity
EToV = [1:K ; 2:K+1]';
return
