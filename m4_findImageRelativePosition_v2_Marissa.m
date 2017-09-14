% zheshigeiwodenvshenxiede. suiranjiandan, danshiwoxiangweitazuoxieshiqing.
%
% determines the relative position of images.  This will generate an input to stitch DIC data.
% Chen Zhe, 2017-04-03
% chenzhe, 2017-04-10 update
% chenzhe, 2017-04-12 update make 'r','c' zero-based for input clarity.


DRAWFIGURE = 0;     % DRAWFIGURE = 0: no figure, 1: some, 2: more
AlignPosition = 1;  % 1 = align mid-point, 2 = align upper left corner
poolSize = 10;      % initial pool size, use to coarsen/average/smooth data to compensate for image distortion in order to do correlation
delta = 100;   % compensate for the distortion on left edge
% path for images and dic data (dic data usually in the same folder as the images)
path_DIC = uigetdir('D:\Marissa_test_20170410_cropped\20170409_ts5Al_01_e0','select the folder that contains the images')
% file name format: [fileNamePrefix_1,iE,fileNamePrefix_2='_', 'r', iR, 'c', iC]
% e.g., 20170409_ts5Al_01_e4_r0c0
fileNamePrefix_1 = '20170409_ts5Al_01_e';
fileNamePrefix_2 = '_';


row_start = 0; % starting # of FOV rows
row_end = 3;
col_start = 0;
col_end = 13;  % ending # of FOV cols

clear transX; clear transY; clear transX_incremental; clear transY_incremental;
transX_incremental(1,1) = 0;
transY_incremental(1,1) = 0;
transX(1,1) = 0;
transY(1,1) = 0;

%%
iE = 5;        % 'e#' in the file name, i.e., stop/pause #  ----------------------------------------------

clear specialRC;
% [iR,iC, iR2,iC2, transX_intremental, transY_incremental]
switch iE
    case 3
        % for e3
        specialRC = [3,1, 3,2,4885,-30;
            3,2, 3,3, 4854,32;
            3,3, 3,4, 4840,-201];  % use this to manually/artificially adjust for the wrong mag figure.
    case 4
        % for e4
        specialRC = [3,6, 3,7, 5874,-316];
    case 5
        % for e5
        specialRC = [0,7, 0,8, 6209,120;
            1,6, 1,7, 6453,-240;
            2,6, 2,7, 6209,0];
end

for iR = row_start:row_end
    for iC = col_start:col_end-1
        pause(1);
        close all;
        fName1 = [fileNamePrefix_1,num2str(iE),fileNamePrefix_2,'r',num2str(iR),'c',num2str(iC)]
        fName2 = [fileNamePrefix_1,num2str(iE),fileNamePrefix_2,'r',num2str(iR),'c',num2str(iC+1)]
        I = imread([path_DIC,'\',fName1,'.tif']);
        J = imread([path_DIC,'\',fName2,'.tif']);
        
        % For debug
        if DRAWFIGURE > 1
            f1 = figure;imshowpair(I,J,'montage');
        end
        
        %         Ipool = my_pool(I,poolSize,1);
        %         Jpool = my_pool(Jprime,poolSize,1);
        
        try
            % initially crop a small region to detect
            xi = floor(size(I,2) * 0.02);
            if iR==0
                yi = floor(size(I,1) * 0.88);
            else
                yi = floor(size(I,1) * 0.02);
            end
            
            xf = xi + floor(size(I,2)*0.1);
            yf = yi + floor(size(I,1)*0.1);
            
            Jprime = J(yi:yf,xi:xf);
            
            c = normxcorr2(Jprime, I);
            if DRAWFIGURE > 1
                figure;surf(c);shading flat; set(gca,'ydir','reverse');view(0,90);
            end
            [ypeak, xpeak] = find(c == max(c(:)));
            xOffSet = xpeak-size(Jprime,2);
            yOffSet = ypeak-size(Jprime,1);
            
            if DRAWFIGURE > 0
                figure;imshowpair(I,J,'montage');
                imrect(gca, [xOffSet+1, yOffSet+1, size(Jprime,2), size(Jprime,1)]);
                imrect(gca, [xi+size(I,2), yi, size(Jprime,2), size(Jprime,1)]);
            end
            
            if (iC>0)&&((xOffSet-(xi-1)<0)||abs(yOffSet-(yi-1))>500)
                error();
            end
            
        catch
            % get rectangular region manually
            try close(f1); catch ; end
            f1 = figure;set(f1,'name','select an area to match');imshowpair(I,J,'montage');
            rect = round(getrect(gcf));
            xi = rect(1); yi = rect(2); xf = xi + rect(3); yf = yi + rect(4);
            xi = mod(xi-size(I,2),size(I,2));
            xf = mod(xf-size(I,2),size(I,2));
            yi = mod(yi-size(I,1),size(I,1));
            yf = mod(yf-size(I,1),size(I,1));
            Jprime = J(yi:yf,xi:xf);
            
            c = normxcorr2(Jprime, I);
            if DRAWFIGURE > 1
                figure;surf(c);shading flat; set(gca,'ydir','reverse');view(0,90);
            end
            [ypeak, xpeak] = find(c == max(c(:)));
            xOffSet = xpeak-size(Jprime,2);
            yOffSet = ypeak-size(Jprime,1);
            
            if DRAWFIGURE > 0
                figure;imshowpair(I,J,'montage');
                imrect(gca, [xOffSet+1, yOffSet+1, size(Jprime,2), size(Jprime,1)]);
                imrect(gca, [xi+size(I,2), yi, size(Jprime,2), size(Jprime,1)]);
            end
            
            c = normxcorr2(Jprime, I);
            
            if DRAWFIGURE > 1
                figure;surf(c);shading flat; set(gca,'ydir','reverse');view(0,90);
            end
            [ypeak, xpeak] = find(c == max(c(:)));
            xOffSet = xpeak-size(Jprime,2);
            yOffSet = ypeak-size(Jprime,1);
            
            if DRAWFIGURE > 0
                figure;imshowpair(I,J,'montage');
                imrect(gca, [xOffSet+1, yOffSet+1, size(Jprime,2)+2*delta, size(Jprime,1)+2*delta]);
            end
            
        end
        
        transX_incremental(iR+1,iC+2) = xOffSet - (xi-1)
        transY_incremental(iR+1,iC+2) = yOffSet - (yi-1)
        transX(iR+1,iC+2) = xOffSet - (xi-1) + transX(iR+1,iC+1);
        transY(iR+1,iC+2) = yOffSet - (yi-1) + transY(iR+1,iC+1);
        
        if (iC==0)&&(iR<row_end)
            pause(1);
            close all;
            fName1 = [fileNamePrefix_1,num2str(iE),fileNamePrefix_2,'r',num2str(iR),'c',num2str(iC)]
            fName3 = [fileNamePrefix_1,num2str(iE),fileNamePrefix_2,'r',num2str(iR+1),'c',num2str(iC)]
            I = imread([path_DIC,'\',fName1,'.tif']);
            J = imread([path_DIC,'\',fName3,'.tif']);
            
            % For debug
            if DRAWFIGURE > 1
                figure;imshowpair(I,J,'montage');
            end
          
            try
                % initially crop a small region to detect
                xi = floor(size(I,2) * 0.45);
                yi = floor(size(I,1) * 0.02);
                xf = xi + floor(size(I,2)*0.1);
                yf = yi + floor(size(I,1)*0.1);
                
                Jprime = J(yi:yf,xi:xf);
                
                c = normxcorr2(Jprime, I);
                if DRAWFIGURE > 1
                    figure;surf(c);shading flat; set(gca,'ydir','reverse');view(0,90);
                end
                [ypeak, xpeak] = find(c == max(c(:)));
                xOffSet = xpeak-size(Jprime,2);
                yOffSet = ypeak-size(Jprime,1);
                
                if DRAWFIGURE > 0
                    figure;imshowpair(I,J,'montage');
                    imrect(gca, [xOffSet+1, yOffSet+1, size(Jprime,2), size(Jprime,1)]);
                    imrect(gca, [xi+size(I,2), yi, size(Jprime,2), size(Jprime,1)]);
                end
                
                if yOffSet-(yi-1)<0
                    error();
                end
                
            catch
                % get rectangular region manually
                try close(f1); catch ; end
                f1 = figure;set(f1,'name','select an area to match');imshowpair(I,J,'montage');
                rect = round(getrect(gcf));
                xi = rect(1); yi = rect(2); xf = xi + rect(3); yf = yi + rect(4);
                xi = mod(xi-size(I,2),size(I,2));
                xf = mod(xf-size(I,2),size(I,2));
                yi = mod(yi-size(I,1),size(I,1));
                yf = mod(yf-size(I,1),size(I,1));
                Jprime = J(yi:yf,xi:xf);
                
                c = normxcorr2(Jprime, I);
                if DRAWFIGURE > 1
                    figure;surf(c);shading flat; set(gca,'ydir','reverse');view(0,90);
                end
                [ypeak, xpeak] = find(c == max(c(:)));
                xOffSet = xpeak-size(Jprime,2);
                yOffSet = ypeak-size(Jprime,1);
                
                if DRAWFIGURE > 0
                    figure;imshowpair(I,J,'montage');
                    imrect(gca, [xOffSet+1, yOffSet+1, size(Jprime,2), size(Jprime,1)]);
                    imrect(gca, [xi+size(I,2), yi, size(Jprime,2), size(Jprime,1)]);
                end
                
                c = normxcorr2(Jprime, I);
                
                if DRAWFIGURE > 1
                    figure;surf(c);shading flat; set(gca,'ydir','reverse');view(0,90);
                end
                [ypeak, xpeak] = find(c == max(c(:)));
                xOffSet = xpeak-size(Jprime,2);
                yOffSet = ypeak-size(Jprime,1);
                
                if DRAWFIGURE > 0
                    figure;imshowpair(I,J,'montage');
                    imrect(gca, [xOffSet+1, yOffSet+1, size(Jprime,2)+2*delta, size(Jprime,1)+2*delta]);
                end
                
            end
            
            transX_incremental(iR+2,iC+1) = xOffSet - (xi-1)
            transY_incremental(iR+2,iC+1) = yOffSet - (yi-1)
            transX(iR+2,iC+1) = xOffSet - (xi-1) + transX(iR+1,iC+1);
            transY(iR+2,iC+1) = yOffSet - (yi-1) + transY(iR+1,iC+1);

        end
        
        % if have specialRC to handle
        if exist('specialRC','var')&&(~isempty(specialRC))
            ind = find((iR==specialRC(:,1))&(iC==specialRC(:,2))&(iR==specialRC(:,3))&(iC+1==specialRC(:,4)));
        else
            ind = [];
        end
        if ~isempty(ind)
            transX_incremental(1+specialRC(ind,3), 1+specialRC(ind,4)) = specialRC(ind,5)
            transY_incremental(1+specialRC(ind,3), 1+specialRC(ind,4)) = specialRC(ind,6)
            transX(1+specialRC(ind,3), 1+specialRC(ind,4)) = specialRC(ind,5) + transX(iR+1,iC+1);
            transY(1+specialRC(ind,3), 1+specialRC(ind,4)) = specialRC(ind,6) + transY(iR+1,iC+1);
        end
        
    end
end



save(['translations_stop_',num2str(iE)],'transX','transY','transX_incremental','transY_incremental');







