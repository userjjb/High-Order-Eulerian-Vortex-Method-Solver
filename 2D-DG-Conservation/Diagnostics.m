clear all
clc

it=1;
tests=30;
for i=tests;
    data{it}=num2str(i);
    it=it+1;
end

post='GDpt4sg_24';
pre= 'K4_';

for runs=1:length(data)
    load([pre,data{runs},post,'.mat'])
    N=setup(2); M=setup(3); del=setup(4); delt=setup(5); K=setup(7:8); B=setup(9:12);
    run('CalcedParams');
    QwPre=(delX/2)^2*reshape(Qw'*Qw,1,[]);
    
    w_elem=reshape(permute(reshape(wxt,Np,K(2),Np,K(1),[]),[1 3 2 4 5]),1,Np^2,K(2)*K(1),[]); %Reshaped to col-wise element chunks
    w_tot=abs(permute(mtimesx(w_elem,QwPre'),[3 4 1 2])); %Sum of vorticity in each elem
    dw(runs,:)=(setup(1)-sum(w_tot,1))/setup(1);
    r_xE=reshape(wxm(Nnumy(:,:,Estreamy)),Np^2,[]);
    r_yE=reshape(wym(Nnumy(:,:,Estreamy)),Np^2,[]);
    Moment=@(mm,nn) reshape(sum(mtimesx(QwPre,bsxfun(@times,permute(w_elem,[2 3 4 1]),r_xE.^mm.*r_yE.^nn)),2),1,[]);
    J20(runs,:)= Moment(2,0); J02(runs,:)= Moment(0,2);
    J11= Moment(1,1); J01(runs,:)= Moment(0,1); J10(runs,:)= Moment(1,0);
    J= J20(runs,:)+J02(runs,:); D= J20(runs,:)-J02(runs,:); R= sqrt(D.^2+4*J11.^2); AspectRatio(runs,:)= sqrt((J+R)./(J-R));
    
    tt(runs,:)= ([1:size(wxt,4)]-1)*delt;
end

rvals= linspace(0,1,length(data));
gvals= linspace(1,0,length(data));
bvals= 0.5*mod((1:length(data))-1,2);
figure(1)
hold on
for runs= 1:length(data)
    hp(runs)= plot(tt(runs,2:end),log10(dw(runs,2:end)),'-k');%'Color',[rvals(runs), gvals(runs),bvals(runs)],'linewidth',2);
end
hl=legend(hp,data,'Location','northwest');
set(hl, 'Interpreter', 'none')

%Plot convergence
figure(2)
hold on
for i=1:size(dw,2)
    C(i,:)=polyfit(log(tests),log(dw(:,i))',1);
end
plot(tt(1,:),C(:,1),'k')
xlabel('Time(s)')
ylabel('Rate of Convergence')
title('Conservation of Vorticity')

figure(3)
hold on
for runs= 1:length(data)
    hp(runs)= plot(tt(runs,:),(J01(runs,:)-J01(runs,1))/J01(runs,1),'-k');%'Color',[rvals(runs), gvals(runs),bvals(runs)],'linewidth',2);
end
hl=legend(hp,data,'Location','northwest');
set(hl, 'Interpreter', 'none')

% figure(4)
% hold on
% for i=1:size(dw,2)
%     CJ(i,:)=polyfit(log(tests),log(J01(:,i)-J01(:,1))',1);
% end
% plot(tt(1,:),CJ(:,1))

% figure(5)
% hold on
% for runs= 1:length(data)
%     hp(runs)= plot(tt(runs,:),AspectRatio(runs,:),'Color',[rvals(runs), gvals(runs),bvals(runs)],'linewidth',2);
% end
% hl=legend(hp,data,'Location','northwest');
% set(hl, 'Interpreter', 'none')