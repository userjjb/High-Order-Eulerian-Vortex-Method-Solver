%Some terminology: the structured tensor grid means we can completely
%separate x and y components and discretize independently. They only couple
%when we calculate dw_p/dt at a point 'p' from the sum of each directions
%contribution, so dw_p/dt = dw_px/dt + dw_py/dt
%The set of colinear points along a coordinate direction will be called a
%"stream"
  
NewEndTime=144.5;
PrevElapsed= setup(end);
PrevEndTime=t+delt;
%Solver--------------------------------------------------------------------
tic
for t=PrevEndTime:delt:NewEndTime
    if mod(StepNum,LogPeriod)==0
        run('PlotNSave')
    end; StepNum= StepNum+1;

    for i=1:nS
        St= t+RKc(i)*delt;              %Unused currently, St is the stage time if needed

%---Velocity eval of current timestep's vorticity config-----------
    v_xB(:)=0; v_yB(:)=0; v_xBF(:)=0; v_yBF(:)=0; v_xE(:)=0; v_yE(:)=0;
    w_elem=reshape(permute(reshape(wy,Np,K(2),Np,K(1)),[1 3 2 4]),1,Np^2,K(2)*K(1)); %Reshaped to col-wise element chunks
    w_tot=abs(permute(mtimesx(w_elem,QwPre'),[3 1 2])); %Sum of vorticity in each elem
    mask=find(w_tot>w_thresh); %Find "important" elements
    w_elemPre=bsxfun(@times,QwPre,w_elem(:,:,mask)); %Pre-multiply by quad weights for speed
    
    v_xI= reshape(mtimesx(w_elemPre, kernel_x),1,Mp-2,Np*(2*NearRange+1)^2,[]);
    v_yI= reshape(mtimesx(w_elemPre, kernel_y),1,Mp-2,Np*(2*NearRange+1)^2,[]);
    for it=1:length(mask)
        w_source=w_elemPre(:,:,it);
        Src= mask(it);
        NsxS=Nsx(1:numS(Src),Src);
        NsyS=Nsy(1:numS(Src),Src);
        %Form specific source kernel by transforming the generalized source
        %kernel to the specific source loc
        kernel_xB= gkernel_xB(:, [1:Np*K(2)] +Np*(K(2)-Enumy(Src)), [1:K(1)+1] +(K(1)-Enumx(Src)) );
        kernel_yB= gkernel_yB(:, [1:Np*K(1)] +Np*(K(1)-Enumx(Src)), [1:K(2)+1] +(K(2)-Enumy(Src)) );
        %Calculate boundary velocities
        v_xBt= permute(mtimesx(w_source,kernel_xB),[2 3 1]); v_xB= v_xB + v_xBt;
        v_yBt= permute(mtimesx(w_source,kernel_yB),[3 2 1]); v_yB= v_yB + v_yBt;
        %Form far field boundary velocities due to source, add to
        %existing far field velocities. Be sure to leave out near-field
        %boundary velocities as these will be included in the whole
        %element evals
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
        
        k2= RKa(i)*k2 + delt*(wx_dt+wy_dt);
        wx= wx+RKb(i)*k2;
        wy= reshape(reshape(wx,K(1)*Np,[])',Np,1,[]); %Reshape wx to match global node ordering
    end
end
setup(end)=toc+PrevElapsed
if saveQ; save(filename,'wxt','setup'); end