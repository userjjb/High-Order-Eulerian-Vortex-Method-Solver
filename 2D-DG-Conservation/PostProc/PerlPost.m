clear all
clc

it=1;
tests=[13];
for i=tests;
    data{it}=num2str(i);
    it=it+1;
end
post='GD';
pre= 'PS24P3_';

for runs=1:length(data)
    load([pre,data{runs},post,'.mat'])
    N=setup(2); M=setup(3); del=setup(4); delt=setup(5); K=setup(7:8); B=setup(9:12);
    run('CalcedParams');
    wxA=permute(reshape(permute(wxt(:,1,Estreamx,:),[1 3 2 4]),Np,Np,K(2),K(1),[]),[2 1 3 4 5]);
    wx0=wxA(:,:,:,:,1);
    
    L2(runs,:)=sqrt(sum(sum(mtimesx(Qw,mtimesx(bsxfun(@minus,wxA,wx0).^2, Qw'/(2*K(1))^2)),3),4));
    tt(runs,:)= ([1:size(wxA,5)]-1)*delt;
end

%Plot L2 error
rvals= linspace(0,1,length(data));
gvals= linspace(1,0,length(data));
bvals= 0.5*mod((1:length(data))-1,2);
figure(1)
hold on
for runs= 1:length(data)
    hp(runs)= plot(tt(runs,1:1:end),log10(L2(runs,1:1:end)),'-.k');%'Color',[rvals(runs), gvals(runs),bvals(runs)],'linewidth',2);
end
hl=legend(hp,data,'Location','northwest');
set(hl, 'Interpreter', 'none')
xlabel('Time(s)')
ylabel('log(L^2 Error)')
title('L^2 Error: Vorticity')

%Plot convergence
figure(2)
hold on
for i=1:size(L2,2)
    C(i,:)=polyfit(log10(tests),-log10(L2(:,i))',1);
end
plot(tt(1,1:2:end),C(1:2:end,1),'-k')
xlabel('Time(s)')
ylabel('Order of Convergence')
title('L^2 Error: Vorticity')