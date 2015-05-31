wxt(:,:,:,itt)=wx;
tt(itt)=t;
%surf(wxm,wym,reshape(wx,Np*K(1),Np*K(2))')
%[-[.831,.696,.563,.43,.3,.168,.032],0.099,.227,.364]
contour(wxm,wym,reshape(wx,Np*K(1),Np*K(2))',PlotInt,'linewidth',1);
axis equal
axis(B)
%axis([B,zmin,zmax])
%Residual calc, used to calc the L^2 norm
%GA1*exp(-((mod((wxm+1)/2-t*cx/2,1)*2-1).^2/Ga1+(mod((wym+1)/2-t*cy/2,1)*2-1).^2/Gb1))
R=w-reshape(wx,Np*K(1),Np*K(2))';
text(B(1),B(4)*0.8,zmax*1.2,['Time: ',num2str(t),char(10),...
    'L^2 norm: ',num2str(sqrt(sum(sum(norm_h.*R.^2)))),char(10),...
    'Mask: ',num2str(length(mask)),char(10),...
    'Done in: ',num2str((EndTime-t)*toc/t),char(10),...
    '\Delta\omega: ',num2str(setup(1)-sum(w_tot))]);
% figure(2)
% contour(wxm,wym,reshape(Cx,Np*K(1),Np*K(2))',0.1:.1:1,'linewidth',1);
% axis equal
% axis(B)
drawnow
itt=itt+1;
if toc>BackupSave
    if saveQ;
        fprintf('Saving at: %i\n',round(toc))
        save(filename,'wxt','setup'); 
    end
    BackupSave=BackupSave+1800;
end
if or(any(wx>max(1,10*zmax)),any(wx<min(-1,10*zmin)))
    beep;pause(0.1); beep;
    fprintf('Solution divergence may have occured, type "return" to continue anyway or "dbstop" to eject...\n')
    keyboard
end
if abs(sum(reshape(wx-wxtCOMP,[],1)))<.01
    fprintf('HIT! t=%d',itt)
end