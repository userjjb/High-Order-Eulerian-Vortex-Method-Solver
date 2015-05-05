fh=figure;
set(gcf,'PaperPosition',[0 0 10 10]);
v=[-[.831,.696,.563,.43,.3,.168,.036],0.099,.227,.364];
tic
endt=size(wxt,4);
for t=1:endt
    set(0, 'CurrentFigure', fh);
    clf reset;
    
    set(gca,'LooseInset',get(gca,'TightInset'))
    contour(wxm,wym,reshape(wxt(:,:,:,t),Np*K(1),Np*K(2))',[-.946:.07:.45],'linewidth',1);
    set(gcf,'color','w');
    axis equal;
    axis(B);
    print(gcf,'-r432','-dtiff',['test',num2str(t)]) %4x108 for a 4x1080 image
    if mod(t,20)<.5
        fprintf('%i ',floor((endt-t)*toc/t))
    end
end