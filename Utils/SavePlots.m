w_tot0=setup(1); N=setup(2); M=setup(3); del=setup(4); delt=setup(5); EndTime=setup(6); K=setup(7:8); B=setup(9:12);
run('CalcedParams')

h=waitbar(0,'Saving plots...');
fh=figure;

%v=[-[.831,.696,.563,.43,.3,.168,.036],0.099,.227,.364];
v=[0.005,0.01,.07:.08:1];
tic

endt=size(wxt,4);
set(gcf,'PaperPosition',[0 0 10 10]);
for t=1:2:endt
    set(0, 'CurrentFigure', fh);
    clf reset;
    
    set(gca,'LooseInset',get(gca,'TightInset'))
    contour(wxm,wym,reshape(wxt(:,:,:,t),Np*K(1),Np*K(2))',v,'linewidth',1);
    set(gcf,'color','w');
    axis equal;
    axis(B);
    print(gcf,'-r432','-dtiff',['test',num2str(t)]) %4x108 for a 4x1080 image
    waitbar(t/endt,h,sprintf('Estimated time: %i',floor((endt-t)*toc/t)))
end

% endt=size(wxG,4);
% set(gcf,'PaperPosition',[0 0 16 9]);
% for t=1:endt
%     set(0, 'CurrentFigure', fh);
%     clf reset;
%     subplot(1,2,1,'Position',[0,0,.5,1])
%     contour(wxm,wym,reshape(wxP(:,:,:,t),Np*K(1),Np*K(2))',[-.946:.07:.45],'linewidth',1);
%     axis equal;
%     axis(B);
%     set(gca,'YTickLabel',[],'XTickLabel',[])
% 
%     subplot(1,2,2,'Position',[0.5,0,.5,1])
%     contour(wxm,wym,reshape(wxG(:,:,:,t),Np*K(1),Np*K(2))',[-.946:.07:.45],'linewidth',1);
%     axis equal;
%     axis(B);
%     set(gca,'YTickLabel',[],'XTickLabel',[])
%     print(gcf,'-r480','-dtiff',['test',num2str(t)]) %4x108 for a 4x1080 image
%     waitbar(t/endt,h,sprintf('Estimated time: %i',floor((endt-t)*toc/t)))
% end
    
close(h)