addChenFunction;
% close all;
% load('D:\_6.mat')

folderName = 'D:\Marissa_test_20170430_stitched_DIC\partial_xyuv_1to38_ref1\';
folderName = uigetdir(folderName,'select folder');
folderName = [folderName,'\'];
for iE = 0:6
    
    fileName = ['_',num2str(iE),'_oly_200'];
    % fileName = ['D:\_',num2str(iE)];
    figureTitlePrefix = ['e',num2str(iE),'-'];
    load([folderName,fileName]);
    
    [f,a,c] = myplot(x,y,sigma); title([figureTitlePrefix,'sigma']);caxis([-1 1]);
    set(f,'position',[100 100 1600 350]);print([folderName,get(a.Title,'String')],'-dtiff');close(f);
    
    [f,a,c] = myplot(x,y,u); title([figureTitlePrefix,'u']); ua=quantile(u(~isnan(u(:))),0.01);ub=quantile(u(~isnan(u(:))),0.99);caxis([ua ub]);
    set(f,'position',[100 100 1600 350]);print([folderName,get(a.Title,'String')],'-dtiff');close(f);
    
    [f,a,c] = myplot(x,y,v); title([figureTitlePrefix,'v']); va=quantile(v(~isnan(v(:))),0.01);vb=quantile(v(~isnan(v(:))),0.99);caxis([va-eps vb+eps]);
    set(f,'position',[100 100 1600 350]);print([folderName,get(a.Title,'String')],'-dtiff');close(f);
    
    [f,a,c] = myplot(x,y,exx); title([figureTitlePrefix,'exx']); exx1=quantile(exx(~isnan(exx(:))),0.2); exx2=quantile(exx(~isnan(exx(:))),0.8); caxis([exx1 exx2]);
    set(f,'position',[100 100 1600 350]);print([folderName,get(a.Title,'String')],'-dtiff');close(f);
    try
        [f,a,c] = myplot(x,y,exx_Lagrange); title([figureTitlePrefix,'exx_Lagrange']); exx1=quantile(exx_Lagrange(:),0.2); exx2=quantile(exx_Lagrange(:),0.8); caxis([exx1 exx2]);
        set(f,'position',[100 100 1600 350]);print([folderName,get(a.Title,'String')],'-dtiff');close(f);
    catch
    end
end

% [f,a,c] = myplot(x,y,sigma); title('sigma');caxis([-1 1]);
% [f,a,c] = myplot(x,y,u); title('u-pixel'); ua=quantile(u(:),0.01);ub=quantile(u(:),0.99);caxis([ua ub]);
% [f,a,c] = myplot(x,y,v); title('v-pixel'); va=quantile(v(:),0.01);vb=quantile(v(:),0.99);caxis([va vb]);
% [f,a,c] = myplot(x,y,exx); title('exx'); exx1=quantile(exx(:),0.2); exx2=quantile(exx(:),0.8); caxis([exx1 exx2]);