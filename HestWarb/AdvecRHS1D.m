function [rhsu] = AdvecRHS1D(u,time, a, Nfp, Nfaces, K, vmapM, vmapP, nx, vmapI, mapI, mapO, rx, Dr, LIFT, Fscale, alpha)
%Pass all the params instead of using global vars, much faster than
%redeclaring globals each run
        
% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

% form field differences at faces
alpha=alpha; %allow overide
du = zeros(Nfp*Nfaces,K); 
du(:) = (u(vmapM)-u(vmapP)).*(a*nx(:)-(1-alpha)*abs(a*nx(:)))/2;

% impose boundary condition at x=0
uin = -sin(a*time);
du (mapI) = (u(vmapI)- uin ).*(a*nx(mapI)-(1-alpha)*abs(a*nx(mapI)))/2;
du (mapO) = 0;

% compute right hand sides of the semi-discrete PDE
rhsu = -a*rx.*(Dr*u) + LIFT*(Fscale.*(du));
return
