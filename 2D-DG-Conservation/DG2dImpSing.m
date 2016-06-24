%Some terminology: the structured tensor grid means we can completely
%separate x and y components and discretize independently. They only couple
%when we calculate dw_p/dt at a point 'p' from the sum of each directions
%contribution, so dw_p/dt = dw_px/dt + dw_py/dt
%The set of colinear points along a coordinate direction will be called a
%"stream"

close all
clear all
clc

tests=16;                   %Test parameter range to iterate over

for yam=1:numel(tests)
    clearvars wxt tt
    
filename=['5_',num2str(tests(yam)),'Gpt3PS2NR1_8.mat'];
saveQ=0;                    %Save time history of state to file
%---Global domain initialization (parameters)------------------------------
B= 3*[-1 1 -1 1];           %left, right, bottom, top
K= [tests(yam) tests(yam)]; %Num elements along x,y
%Solver parameters
delt= 1;                  %Timestep
del=.3*((B(2)-B(1))/K(1)); %Unused for BS kernel
N= 5;                       %Local vorticity poly order
M= 5;                       %Local velocity poly order
[RKa,RKb,RKc,nS]= LSRKcoeffs('NRK14C');
w_thresh=3*(48^2/prod(K))*1E-9;
ThresholdMap= 1;            %Plot thresholding of elements
EndTime=194.5;
LogPeriod= uint64(1);
BCtype= 'NoInflow';
KernelType='PS2';
ModKernelFile='Iwm55.mat';  %File containing modified kernel values
NearRange=ceil(K(1)*(3/24)); %6th: 2/18<good enough<3/24, 10/18-full coupling (may be prob/order depend)
TestCases=12:13;
alpha= 1;                   %Numerical flux param (1 upwind,0 CD)
PlotInt=[0.001,0.005,0.01,.07:.08:1];%[10.^[-15:2:-1],0.02:.1:1,.98];

%Calculate all derived solver parameters (node/boundary/element positions
%and numbering, discrete norm, and pre-allocate vorticity/velocity vars
run('CalcedParams')
%Setup initial conditions at t_0
[w, prior_toc]=InitialConditions(w,TestCases,wxm,wym,'none');

%Solver--------------------------------------------------------------------
run('SolverSetup')
tic
for t=0:delt:EndTime %timestep
    if mod(StepNum,LogPeriod)==0
        run('PlotNSave') %Plots current state; autosaves every so often
    end; StepNum= StepNum+1;

    for i=1:nS %stagestep
        St= t+RKc(i)*delt; %Unused currently, St is the stage time if needed
        %%A "short" explanation of velocity related calculations:
        %- We calculate the velocity field in two parts: near and far-field. Only elements containing
        %total vorticity greater than 'w_thresh' are included as sources.
        %- Elements within 'NearRange' of a source element are considered near-field; velocities due
        %to the effect of the vorticity within the source element are calculated for all 
        %interpolation nodes.
        %- Elements outside of 'NearRange' from the source element are considered far-field, and the
        %velocities due to the source element are only calculated on the boundaries of far-field 
        %elements.
        %- The nodal stiffness eval is then the sum of the near and far-field evals, at high and low
        %fidelity respectively, using the modified quadrature matrices 'QwSM' and 'QwSMlow'
        %respectively.
        %- To avoid excessive memory indexing we calculate ALL boundary velocities (on x-lines for 
        %example) due to a source: 'v_xBt', the portion of 'v_xBt' outside of 'NearRange' of the 
        %source element forms the far-field boundary velocities due to that source: v_xBFt. If we 
        %loop over all sources and add all 'v_xBt' together we get all of the TOTAL boundary
        %velocities: 'v_xB'. These are used to calculate the boundary fluxes. If we loop over all 
        %sources and add all 'v_xBFt' together we get the total boundary velocities due to ONLY
        %far-field effects: 'v_xBF'.
        %- We also need to calculate the high-fidelity total velocities due to near-field effects.
        %The near field is typically much smaller than the whole domain, so the kernel values for
        %all elements near all sources can be precalculated and stored. This needs only one 
        %evaluation of near-field internal elemental velocities for ALL source elements, this is
        %much quicker than the per-source looping required for the boundary velocities.
        %- To avoid recalulating boundary velocities that have been previously calculated only the 
        %internal near-field velocities are calculated: 'v_xI'. The total element velocities are 
        %reconstructed by combining the near-field boundary velocities from 'v_xBt' and the 
        %near-field internal velocities 'v_xI' inside of the per-source loop to form the elemental 
        %near-field velocities: 'v_xE'.
        %- Ultimately 'v_xB' is used to calculate boundary fluxes, 'v_xBF' to calculate far-field
        %nodal stiffness, and 'v_xE' to calculate near-field nodal stiffness.
        %% ---Velocity eval of current timestep's vorticity config-----------
        v_xB(:)=0; v_yB(:)=0; v_xBF(:)=0; v_yBF(:)=0; v_xE(:)=0; v_yE(:)=0;
        w_elem= reshape(permute(reshape(wy,Np,K(2),Np,K(1)),[1 3 2 4]),1,Np^2,K(2)*K(1)); %Reshaped to col-wise element chunks
        w_tot= abs(permute(mtimesx(w_elem,QwPre),[3 1 2]));%Sum of vorticity in each elem
        mask= find(w_tot>w_thresh);                         %Find "important" elements
        w_elemPre= bsxfun(@times,QwPre',w_elem(:,:,mask));   %Pre-multiply by quad weights for speed
        %Internal element only near-field velocities for ALL sources
        v_xI= reshape(mtimesx(w_elemPre, kernel_x),1,Mp-2,Np*(2*NearRange+1)^2,[]);
        v_yI= reshape(mtimesx(w_elemPre, kernel_y),1,Mp-2,Np*(2*NearRange+1)^2,[]);
        for it=1:length(mask)
            w_source=w_elemPre(:,:,it);
            Src= mask(it);
            %Near streams for this source
            NsxS=Nsx(1:numS(Src),Src);
            NsyS=Nsy(1:numS(Src),Src);
            %Form specific source kernel by transforming the generalized source
            %kernel to the specific source loc
            kernel_xB= gkernel_xB(:, [1:Np*K(2)] +Np*(K(2)-Enumy(Src)), [1:K(1)+1] +(K(1)-Enumx(Src)) );
            kernel_yB= gkernel_yB(:, [1:Np*K(1)] +Np*(K(1)-Enumx(Src)), [1:K(2)+1] +(K(2)-Enumy(Src)) );
            %Calculate boundary velocities
            v_xBt= permute(mtimesx(w_source,kernel_xB),[2 3 1]); v_xB= v_xB + v_xBt;
            v_yBt= permute(mtimesx(w_source,kernel_yB),[3 2 1]); v_yB= v_yB + v_yBt;
            %Form far field boundary velocities due to source, add to existing far field velocities.
            %Be sure to leave out near-field boundary velocities as these will be included in the 
            %whole element evals
            v_xBFt= [v_xBt(EBl),v_xBt(EBr)]; v_xBFt(1,:,NsxS)=0; v_xBF= v_xBF + v_xBFt;
            v_yBFt= [v_yBt(EBb),v_yBt(EBt)]; v_yBFt(1,:,NsyS)=0; v_yBF= v_yBF + v_yBFt;
            %Assemble elementwise velocities for elements nearby the source
            v_xE(1,:,NsxS)= v_xE(1,:,NsxS)+ [v_xBt(EBl(NsxS)), v_xI(1,:,Lsx(1:numS(Src),Src),it) ,v_xBt(EBr(NsxS))];
            v_yE(1,:,NsyS)= v_yE(1,:,NsyS)+ [v_yBt(EBb(NsyS)), v_yI(1,:,Lsy(1:numS(Src),Src),it) ,v_yBt(EBt(NsyS))];
        end
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
        
        k2= RKa(i)*k2 + delt*(wx_dt+wy_dt); %No need to reset k2 at stage start, as RKa(1)=0
        wx= wx+RKb(i)*k2;
        wy= reshape(reshape(wx,K(1)*Np,[])',Np,1,[]); %Reshape wx to match global node ordering
    end
end
setup(16)=toc+prior_toc;
if saveQ; save(filename,'wxt','setup'); end
end