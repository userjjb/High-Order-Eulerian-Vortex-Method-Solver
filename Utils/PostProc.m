clear all
close all
clc

load('6_48_100_Stage.mat')
wxG=wxt;
load('6_48_100_Step.mat')
wxP=wxt;

clearvars wxt
w_tot0=setup(1); N=setup(2); M=setup(3); del=setup(4); delt=setup(5); EndTime=setup(6); K=setup(7:8); B=setup(9:12);
run('CalcedParams')
QwPre=(delX/2)^2*reshape(Qw'*Qw,1,[]);

w_elem= reshape(permute(reshape(wxG,Np,K(2),Np,K(1),[]),[1 3 2 4 5]),1,Np^2,K(2)*K(1),[]);
w_totG= sum(permute(mtimesx(w_elem,QwPre'),[3 4 1 2]));
w_elem= reshape(permute(reshape(wxP,Np,K(2),Np,K(1),[]),[1 3 2 4 5]),1,Np^2,K(2)*K(1),[]);
w_totP= sum(permute(mtimesx(w_elem,QwPre'),[3 4 1 2]));

w_elem= reshape(permute(reshape((wxG-wxP).^2,Np,K(2),Np,K(1),[]),[1 3 2 4 5]),1,Np^2,K(2)*K(1),[]);
w_totL2= sqrt(sum(permute(mtimesx(w_elem,QwPre'),[3 4 1 2])));
clearvars w_elem

plot(w_totG,'g')
hold on
plot(w_totP,'r')

figure(2)
% L2=squeeze(sqrt(sum(sum(bsxfun(@times,norm_h',abs(reshape(wxP,Np*K(1),Np*K(2),[])-reshape(wxG,Np*K(1),Np*K(2),[])).^2)))));
% plot(L2,'r')
% hold on
plot(w_totL2,'b')

%Find G step with closest L2 norm
% w_elem= reshape(permute(reshape(bsxfun(@minus,wxG,wxP(:,:,:,250)).^2,Np,K(2),Np,K(1),[]),[1 3 2 4 5]),1,Np^2,K(2)*K(1),[]);
% w_totL2= sqrt(sum(squeeze(abs(permute(mtimesx(w_elem,QwPre'),[3 1 2 4]))))); plot(w_totL2,'g')

w_elem= reshape(permute(reshape(bsxfun(@minus,wxR,wxP(:,:,:,301)).^2,Np,K(2),Np,K(1),[]),[1 3 2 4 5]),1,Np^2,K(2)*K(1),[]);
w_totL2= sqrt(sum(squeeze(abs(permute(mtimesx(w_elem,QwPre'),[3 1 2 4]))))); plot(w_totL2,'g')
w_elem= reshape(permute(reshape(bsxfun(@minus,wxR,wxG(:,:,:,309)).^2,Np,K(2),Np,K(1),[]),[1 3 2 4 5]),1,Np^2,K(2)*K(1),[]);
w_totL2= sqrt(sum(squeeze(abs(permute(mtimesx(w_elem,QwPre'),[3 1 2 4]))))); plot(w_totL2,'g')