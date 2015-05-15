%Some terminology: the structured tensor grid means we can completely
%separate x and y components and discretize independently. They only couple
%when we calculate dw_p/dt at a point 'p' from the sum of each directions
%contribution, so dw_p/dt = dw_px/dt + dw_py/dt
%The set of colinear points along a coordinate direction will be called a
%"stream"

close all
clear all
clc

%Should match G309 best, compare against P301
filename='6_12_100_Stage_NR5.mat';
saveQ=1;
%Solver parameters
alpha= 1;                           %Numerical flux param (1 upwind,0 CD)
N= 6;                               %Local vorticity poly order
M= 6;                               %Local velocity poly order
[RKa,RKb,RKc,nS]= LSRKcoeffs('NRK14C');
w_thresh=4E-9;
del=0.5*(7.875/20);
delt= 1*.32;
EndTime=100;
LogPeriod= uint64(1);
BCtype= 'NoInflow';
NearRange=5;
TestCases=5:8;
%---Global domain initialization (parameters)------------------------------
B= 3.5*[-1.25 1 -1.25 1];           %left, right, bottom, top
K= [12 12];               %Num elements along x,y

%Calculate all derived solver parameters (node/boundary/element positions
%and numbering, discrete norm, and pre-allocate vorticity/velocity vars
run('CalcedParams2')
%Setup initial conditions at t_0
w=InitialConditions(w,TestCases,wxm,wym);

%Solver--------------------------------------------------------------------
run('SolverSetup2')
tic
for t=0:delt:EndTime
    if mod(StepNum,LogPeriod)==0
        run('PlotNSave')
    end; StepNum= StepNum+1;
    
    for i=1:nS
        St= t+RKc(i)*delt;              %Unused currently, St is the stage time if needed
        
        %---Velocity eval of current timestep's vorticity config-----------
        v_xE(:)=0; v_yE(:)=0;
        w_elem=reshape(permute(reshape(wy,Np,K(2),Np,K(1)),[1 3 2 4]),1,Np^2,K(2)*K(1)); %Reshaped to col-wise element chunks
        w_tot=abs(permute(mtimesx(w_elem,QwPre'),[3 1 2])); %Sum of vorticity in each elem
        mask=find(w_tot>w_thresh); %Find "important" elements
        w_elemPre=bsxfun(@times,QwPre,w_elem); %Pre-multiply by quad weights for speed

        v_xI= reshape(mtimesx(w_elemPre(:,:,mask), kernel_x),1,Mp-1,Np*(2*NearRange+1)^2,[]);
        v_yI= reshape(mtimesx(w_elemPre(:,:,mask), kernel_y),1,Mp-1,Np*(2*NearRange+1)^2,[]);
        %Calculate boundary velocities
        v_xB= mtimesx(w_elemPre(:),'T',kernel_xB);
        v_yB= mtimesx(w_elemPre(:),'T',kernel_yB);
        v_xBF= [v_xB(EBl),v_xB(EBr)];
        v_yBF= [v_yB(EBb),v_yB(EBt)];
        for it=1:length(mask)
            Src= mask(it);
            %Assemble elementwise velocities for elements nearby the source
            v_xE(1,:,Nsx(1:numS(Src),Src))=...
                v_xE(1,:,Nsx(1:numS(Src),Src)) + v_xI(1,:,Lsx(1:numS(Src),Src),it);
            v_yE(1,:,Nsy(1:numS(Src),Src))=...
                v_yE(1,:,Nsy(1:numS(Src),Src)) + v_yI(1,:,Lsy(1:numS(Src),Src),it);
        end
        v_xE(1,end+1,:)=%!!!!!!TODO: Add right boundary velocity
        v_yE(1,end+1,:)=%^^^^^^^^^^^^^^^^top
        
        %!!!!!!!TODO: Add in right/top global domain boundary velocities absent in
        %v_xE(EBr) / v_yE(EBt)
        v_xB= reshape( v_xB' + squeeze(v_xE(1,1,:)) ,Np*K(2),K(1)+1 );
        v_yB= reshape( v_yB' + squeeze(v_yE(1,1,:)) ,K(2)+1,Np*K(1) );
        %---Velocity eval ends---------------------------------------------
        
        %---Advection------------------------------------------------------
        w_lx= mtimesx(Ll',wx);          %Left interpolated vorticity
        w_rx= mtimesx(Lr',wx);          %Right interpolated vorticity
        w_bx= mtimesx(Ll',wy);          %Bottom interpolated vorticity
        w_tx= mtimesx(Lr',wy);          %Top interpolated vorticity
        
        %Boundary fluxes
        if BCtype== 'NoInflow'
            v_xBC= v_xB; v_xBC(:,1)= min(v_xBC(:,1),0); v_xBC(:,end)= max(v_xBC(:,end),0);
            v_yBC= v_yB; v_yBC(1,:)= min(v_yBC(1,:),0); v_yBC(end,:)= max(v_yBC(end,:),0);
        end
        fl= abs( v_xBC(EBl) ).*( w_rx(x_km1).*(sign(v_xB(EBl))+alpha) + w_lx.*(sign(v_xB(EBl))-alpha) );
        fr= abs( v_xBC(EBr) ).*( w_rx.*(sign(v_xB(EBr))+alpha) + w_lx(x_kp1).*(sign(v_xB(EBr))-alpha) );
        fb= abs( v_yBC(EBb) ).*( w_tx(y_km1).*(sign(v_yB(EBb))+alpha) + w_bx.*(sign(v_yB(EBb))-alpha) );
        ft= abs( v_yBC(EBt) ).*( w_tx.*(sign(v_yB(EBt))+alpha) + w_bx(y_kp1).*(sign(v_yB(EBt))-alpha) );
        
        %Nodal total surface flux
        SurfFlux_x=bsxfun(@times,fr,LrM)-bsxfun(@times,fl,LlM);
        SurfFlux_y=bsxfun(@times,ft,LrM)-bsxfun(@times,fb,LlM);
        %Nodal stiffness eval
        Stiff_x= mtimesx(v_xBF,mtimesx(QwSMlow,wx));
        Stiff_x= Stiff_x + mtimesx(v_xE,mtimesx(QwSM,wx));
        Stiff_y= mtimesx(v_yBF,mtimesx(QwSMlow,wy));
        Stiff_y= Stiff_y + mtimesx(v_yE,mtimesx(QwSM,wy));

        wx_dt= permute(Stiff_x-SurfFlux_x,[4 1 3 2]); %Reshape to match wx
        wy_dt= reshape(reshape(Stiff_y-SurfFlux_y,K(2),[])',Np,1,[]); %Reshape to match wx
        
        k2= RKa(i)*k2 + delt*(wx_dt+wy_dt);
        wx= wx+RKb(i)*k2;
        wy= reshape(reshape(wx,K(1)*Np,[])',Np,1,[]); %Reshape wx to match global node ordering
    end
end
if saveQ; save(filename,'wxt','setup'); end