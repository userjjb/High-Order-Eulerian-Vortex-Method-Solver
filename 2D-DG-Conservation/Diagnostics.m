clear all
clc

it=1;
tests=[3:2:13];
for i=tests;
    data{it}=num2str(i);
    it=it+1;
end
post='GDpt4ps2';
pre= 'P3_';

for runs=1:length(data)
    load([pre,data{runs},post,'.mat'])
    N=setup(2); M=setup(3); del=setup(4); delt=setup(5); K=setup(7:8); B=setup(9:12);
    run('CalcedParams');
    QwPre=(delX/2)^2*reshape(Qw'*Qw,1,[]);
    
    w_elem=reshape(permute(reshape(wxt,Np,K(2),Np,K(1),[]),[1 3 2 4 5]),1,Np^2,K(2)*K(1),[]); %Reshaped to col-wise element chunks
    w_tot=abs(permute(mtimesx(w_elem,QwPre'),[3 4 1 2])); %Sum of vorticity in each elem
    dw(runs,:)=(setup(1)-sum(w_tot))/setup(1);
    
    tt(runs,:)= ([1:size(wxt,4)]-1)*delt;
end

%Plot L2 error
rvals= linspace(0,1,length(data));
gvals= linspace(1,0,length(data));
bvals= 0.5*mod((1:length(data))-1,2);
figure(1)
hold on
for runs= 1:length(data)
    hp(runs)= plot(tt(runs,:),log10(dw(runs,:)),'Color',[rvals(runs), gvals(runs),bvals(runs)],'linewidth',2);
end
hl=legend(hp,data,'Location','northwest');
set(hl, 'Interpreter', 'none')

%Plot convergence
figure(2)
hold on
for i=1:size(dw,2)
    C(i,:)=polyfit(log(tests),log(dw(:,i))',1);
end
plot(tt(1,:),C(:,1))