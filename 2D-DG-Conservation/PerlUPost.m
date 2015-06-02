clear all
clc

vAfx=@(x,y) y.*( (-1./(16*(x.^2+y.^2))) .* (1-(1-min(x.^2+y.^2,1)).^8) );

it=1;
tests=[1:7];
for i=tests;
    data{it}=num2str(i);
    it=it+1;
end
post='GD';
pre= 'PS1P6_';

for runs=1:length(data)
    clearvars v_xEE
    load([pre,data{runs},post,'.mat'])
    N=setup(2); M=setup(3); del=setup(4); delt=setup(5); K=setup(7:8); B=setup(9:12);
    TestCases=0;
    EndTime=0;
    w_thresh=0;
    BCtype= 'NoInflow';
    KernelType='PS2';
    NearRange=ceil(K(1)-1);
    
    run('CalcedParams');
    run('SolverSetup')
    
    x_vE=reshape(bsxfun(@plus,(Qx2+1)*(delX/2),Ex(1:end-1)),1,[]);
    [t2,t1]= meshgrid(y_w,x_vE);
    rv_xE= reshape([t1(:),t2(:)]',1,2,Mp,[]);
    
    v_xA= vAfx(rv_xE(1,1,:,:), rv_xE(1,2,:,:));
    v_xA(isnan(v_xA)) = 0; %Fix any NaNs occuring at the singularity
    v_xAA= permute(reshape(permute(v_xA(1,1,:,Estreamx),[3 4 1 2]),Np,Mp,K(2),K(1)),[2 1 3 4]);
    
    for t=1:size(wxt,4)
    v_xB(:)=0; v_xBF(:)=0; v_xE(:)=0;
    w_elem=reshape(permute(reshape(wxt(:,:,:,t),Np,K(2),Np,K(1)),[1 3 2 4]),1,Np^2,K(2)*K(1)); %Reshaped to col-wise element chunks
    w_tot=abs(permute(mtimesx(w_elem,QwPre'),[3 1 2])); %Sum of vorticity in each elem
    mask=find(w_tot>w_thresh); %Find "important" elements
    w_elemPre=bsxfun(@times,QwPre,w_elem(:,:,mask)); %Pre-multiply by quad weights for speed
    
    v_xI= reshape(mtimesx(w_elemPre, kernel_x),1,Mp-2,Np*(2*NearRange+1)^2,[]);
    for it=1:length(mask)
        w_source=w_elemPre(:,:,it);
        Src= mask(it);
        NsxS=Nsx(1:numS(Src),Src);
        %Form specific source kernel by transforming the generalized source
        %kernel to the specific source loc
        kernel_xB= gkernel_xB(:, [1:Np*K(2)] +Np*(K(2)-Enumy(Src)), [1:K(1)+1] +(K(1)-Enumx(Src)) );
        %Calculate boundary velocities
        v_xBt= permute(mtimesx(w_source,kernel_xB),[2 3 1]); v_xB= v_xB + v_xBt;
        %Form far field boundary velocities due to source, add to
        %existing far field velocities. Be sure to leave out near-field
        %boundary velocities as these will be included in the whole
        %element evals
        v_xBFt= [v_xBt(EBl),v_xBt(EBr)]; v_xBFt(1,:,NsxS)=0; v_xBF= v_xBF + v_xBFt;
        %Assemble elementwise velocities for elements nearby the source
        v_xE(1,:,NsxS)= v_xE(1,:,NsxS)+ [v_xBt(EBl(NsxS)), v_xI(1,:,Lsx(1:numS(Src),Src),it) ,v_xBt(EBr(NsxS))];
    end
    v_xEE(:,:,:,:,t)= permute(reshape(permute(v_xE(1,:,Estreamx),[2 3 1]),Np,Mp,K(2),K(1)),[2 1 3 4]);    
    end
    L2(runs,:)=sqrt(sum(sum(mtimesx(Qw2,mtimesx(bsxfun(@minus,v_xAA,v_xEE).^2, Qw'/(2*K(1))^2)),3),4));
    tt(runs,:)= ([1:size(wxt,4)]-1)*delt;
end

%Plot L2 error
rvals= linspace(0,1,length(data));
gvals= linspace(1,0,length(data));
bvals= 0.5*mod((1:length(data))-1,2);
figure(1)
hold on
for runs= 1:length(data)
    hp(runs)= plot(tt(runs,:),log10(L2(runs,:)),'-','Color',[rvals(runs), gvals(runs),bvals(runs)],'linewidth',1);
end
hl=legend(hp,data,'Location','northwest');
set(hl, 'Interpreter', 'none')

%Plot convergence
figure(2)
hold on
for i=1:size(L2,2)
    C(i,:)=polyfit(log10(tests),log10(L2(:,i))',1);
end
plot(tt(1,1:2:end),C(1:2:end,1),'-ko')
xlabel('Time(s)')
ylabel('Order of Convergence')
title('L^2 Error: Velocity')