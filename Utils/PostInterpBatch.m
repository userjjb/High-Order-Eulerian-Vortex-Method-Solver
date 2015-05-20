clear all
clc

load('9_36GD.mat')
N=setup(2); M=setup(3); del=setup(4); deltE=setup(5); K=setup(7:8); B=setup(9:12);
run('CalcedParams'); NpE=Np; QxE=Qx; KE=K(1); delXE=delX;
wxE=permute(reshape(permute(wxt(:,1,Estreamx,:),[1 3 2 4]),Np,Np,K(2),K(1),[]),[2 1 3 4 5]);

ndE=(QxE+1)/2;
nnE=elim(ndE',ndE',[1 3 2]);
LagE= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,nnE(nv,:,:)),bsxfun(@minus,ndE(nv),nnE(nv,:,:))),3);

interp_g=@(x,y,G) LagE(x,1:Np)'*(G*LagE(y,1:Np));
%--------

load('6-2_12GD.mat')
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

freqE= lcm(deltA*100,deltE*100)/(deltE*100);
freqA= lcm(deltA*100,deltE*100)/(deltA*100);
tratio= freqA/freqE;
endtE=size(wxE,5);
endtA=size(wxA,5);

h=waitbar(0);
it=1;
tic
for t=1:freqE:endtE
    wxIA= permute(mtimesx(LxA, mtimesx( wxA(:,:,:,:,1+(t-1)*tratio),LyA )),[1 3 2 4]);
    wxIE= permute(mtimesx(LxE, mtimesx( wxE(:,:,:,:,t),LyE )),[1 3 2 4]);

    L2(it)= sqrt((delXE/IdE)^2*sum((wxIE(:)-wxIA(:)).^2));
    tt(it)= deltE*(t-1);
    waitbar(t/endtE,h,sprintf('Estimated time: %i',floor((endtE-t)*toc/t)))
    it= it+1;
end
close(h)
plot(tt,L2,'r')
hold on

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