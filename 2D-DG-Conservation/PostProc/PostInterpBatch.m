clear all
clc
set(0, 'DefaulttextInterpreter', 'none');

it=1;
% tests=[28];
% for i=tests;
%     data{it}=num2str(i);
%     it=it+1;
% end
data = {'6_18Gpt3PS2NR10' '6_36Gpt3PS2'};
post='';
pre= '';

load('3_108Gpt3PS2NR1_8.mat')
N=setup(2); M=setup(3); del=setup(4); deltE=setup(5); K=setup(7:8); B=setup(9:12);
run('CalcedParams'); NpE=Np; QxE=Qx; KE=K(1); delXE=delX;
wxE=permute(reshape(permute(wxt(:,1,Estreamx,:),[1 3 2 4]),Np,Np,K(2),K(1),[]),[2 1 3 4 5]);

ndE=(QxE+1)/2;
nnE=elim(ndE',ndE',[1 3 2]);
LagE= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,nnE(nv,:,:)),bsxfun(@minus,ndE(nv),nnE(nv,:,:))),3);

interp_g=@(x,y,G) LagE(x,1:Np)'*(G*LagE(y,1:Np));
%--------

for runs=1:length(data)
load([pre,data{runs},post,'.mat'])
N=setup(2); M=setup(3); del=setup(4); deltA=setup(5); K=setup(7:8); B=setup(9:12);
run('CalcedParams'); NpA=Np; QxA=Qx; KA=K(1); delXA=delX;
wxA=permute(reshape(permute(wxt(:,1,Estreamx,:),[1 3 2 4]),Np,Np,K(2),K(1),[]),[2 1 3 4 5]);

ndA=(QxA+1)/2;
nnA=elim(ndA',ndA',[1 3 2]);
LagA= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,nnA(nv,:,:)),bsxfun(@minus,ndA(nv),nnA(nv,:,:))),3);

if (lcm(KA,KE)/KA) > (lcm(KA,KE)/KE)
    IdE= max( lcm(KA,KE)/KE , ceil((NpE*2)/(lcm(KA,KE)/KE))*(lcm(KA,KE)/KE) );
    IdA= (IdE/(lcm(KA,KE)/KE))* (lcm(KA,KE)/KA);
else
    IdA= max( lcm(KA,KE)/KA , ceil((NpA*2)/(lcm(KA,KE)/KA))*(lcm(KA,KE)/KA) );
    IdE= (IdA/(lcm(KA,KE)/KA))* (lcm(KA,KE)/KE);
end

IpE= 1/(2*IdE) : 1/IdE : 1-1/(2*IdE);
IpA= 1/(2*IdA) : 1/IdA : 1-1/(2*IdA);

LxE= LagE(IpE,1:NpE)';
LyE= LagE(IpE,1:NpE);
LxA= LagA(IpA,1:NpA)';
LyA= LagA(IpA,1:NpA);

%Match the timestepping frequencies
freqE= lcm(deltA*1000,deltE*1000)/(deltE*1000);
freqA= lcm(deltA*1000,deltE*1000)/(deltA*1000);
tratio= freqA/freqE;
endtE=size(wxE,5);
endtA=size(wxA,5);

h=waitbar(0);
it=1;
tic
%wxIA= permute(mtimesx(LxA, mtimesx( wxA(:,:,:,:,end-8),LyA )),[1 3 2 4]);
for t=1:freqE:endtE
    tA=1+(t-1)*tratio; if tA>endtA break; end
    wxIA= permute(mtimesx(LxA, mtimesx( wxA(:,:,:,:,tA),LyA )),[1 3 2 4]);
    wxIE= permute(mtimesx(LxE, mtimesx( wxE(:,:,:,:,t),LyE )),[1 3 2 4]);

    L2(runs,it)= sqrt(((delXE/(B(1)-B(2)))/IdE)^2*sum((wxIE(:)-wxIA(:)).^2));
    L2d(runs,it)= sqrt(((delXE/(B(1)-B(2)))/IdE)^2*sum(((wxIE(:)-wxIA(:))/20).^2));
    %L1(runs,it)= (((delXE/(B(1)-B(2)))/IdE)^2*sum(abs(wxIE(:)-wxIA(:))));
    tt(runs,it)= deltE*(t-1);
    waitbar(t/endtE,h,sprintf('%s: %i',[pre,data{runs},post],floor((endtE-t)*toc/t)))
    it= it+1;
end
PlotPoints(runs)= it-1;

end
close(h)

rvals= linspace(0,1,length(data));
gvals= linspace(1,0,length(data));
bvals= 0.5*mod((1:length(data))-1,2);
figure(1)
hold on
for runs= 1:length(data)
    hp(runs)= plot(tt(runs,1:1:PlotPoints(runs)),(L2(runs,1:1:PlotPoints(runs))),'Color',[rvals(runs), gvals(runs),bvals(runs)]);
end
hl=legend(hp,data,'Location','northwest');
set(hl, 'Interpreter', 'tex')

set(0, 'DefaulttextInterpreter', 'tex');

figure(2)
hold on
for i=1:size(L2,2)
    C2(i,:)=polyfit(log(tests),log(L2(:,i))',1);
    %C1(i,:)=polyfit(log(tests),log(L1(:,i))',1);
end
plot(tt(1,:),C2(:,1),'--k')
%plot(tt(1,:),C1(:,1),'r')

% endt=1+round((size(wxE,5)-1)/5);
% tic
% for t=1:endt
%     wxIE= permute(mtimesx(Lx,mtimesx( wxE(:,:,:,:,5*(t-1)+1),Ly )),[1 3 2 4]);
%     wxIA= permute(mtimesx(Lx2, mtimesx( wxA(:,:,:,:,t),Ly2 )),[1 3 2 4]);
% 
%     L2(t)= sqrt((delX/Id)^2*sum((wxIE(:)-wxIA(:)).^2));
%     waitbar(t/endt,h,sprintf('Estimated time: %i',floor((endt-t)*toc/t)))
% end

% IpD=bsxfun(@plus,Ip'*delX,Ex(1:end-1));
% [xx,yy]=meshgrid(IpD(:),IpD(:));
% 
% contour(xx,yy,wxI,[-.946:.07:.45],'linewidth',.3);
% axis equal; axis(B);