% zheshigeiwodenvshenxiede. suiranjiandan, danshiwoxiangweitazuoxieshiqing.
%
% determines the relative position of images.  This will generate an input to stitch DIC data.
% Chen Zhe, 2017-04-03
% chenzhe, 2017-04-10 update


DRAWFIGURE = 0;     % DRAWFIGURE = 0: no figure, 1: some, 2: more
AlignPosition = 1;  % 1 = align mid-point, 2 = align upper left corner
poolSize = 10;      % initial pool size, use to coarsen/average/smooth data to compensate for image distortion in order to do correlation
delta = 100;   % compensate for the distortion on left edge
% path for images and dic data (dic data usually in the same folder as the images)
path_DIC = uigetdir('D:\Marissa_test_20170410_cropped\20170409_ts5Al_01_e0','select the folder that contains the images')
% file name format: [f1,STOP{#},'_',FOV{#,#}]
f1 = '20170409_ts5Al_01_e';
stop_start = 0;     % starting # of stops
stop_end = 10;
f2 = '_';
nFov_row_start = 0; % starting # of FOV rows
nFov_row_end = 3;
nFov_col_start = 0;
nFov_col_end = 13;  % ending # of FOV cols

iStop = 0;
for ii = stop_start:stop_end
    iStop = iStop + 1;
    STOP{iStop} = num2str(ii);
end
iR = 0;
for ii = nFov_row_start:nFov_row_end
    iR = iR + 1;
    iC = 0;
    for jj = nFov_col_start:nFov_col_end
        iC = iC + 1;
        FOV{iR,iC} = ['r',num2str(ii),'c',num2str(jj)];
    end
end
% the 'ugly' way
% STOP = {'0','1','2','3','4','5','6','7','8','9','10'};
% FOV = {
% 'r0c0', 'r0c1', 'r0c2', 'r0c3', 'r0c4', 'r0c5', 'r0c6', 'r0c7', 'r0c8', 'r0c9', 'r0c10', 'r0c11','r0c12', 'r0c13';
% 'r1c0', 'r1c1', 'r1c2', 'r1c3', 'r1c4', 'r1c5', 'r1c6', 'r1c7', 'r1c8', 'r1c9', 'r1c10', 'r1c11','r1c12', 'r1c13';
% 'r2c0', 'r2c1', 'r2c2', 'r2c3', 'r2c4', 'r2c5', 'r2c6', 'r2c7', 'r2c8', 'r2c9', 'r2c10', 'r2c11','r2c12', 'r2c13';
% 'r3c0', 'r3c1', 'r3c2', 'r3c3', 'r3c4', 'r3c5', 'r3c6', 'r3c7', 'r3c8', 'r3c9', 'r3c10', 'r3c11','r3c12', 'r3c13';
% };
clear transX; clear transY; clear transX_incremental; clear transY_incremental;

%%
eNumber = 5;        % 'e#' in the file name, i.e., stop/pause #  ----------------------------------------------

clear specialRC;
% note this is iR,iC. i.e., start from 1. !!!
% [iR,iC, iR2,iC2, transX_intremental, transY_incremental]
switch eNumber
    case 3
% for e3
specialRC = [4,2,4,3,4885,-30;
    4,3, 4,4, 4854,32;
    4,4, 4,5, 4840,-201];  % use this to manually/artificially adjust for the wrong mag figure.
    case 4
% for e4
specialRC = [4,7, 4,8, 5874,-316];
    case 5
% for e5
specialRC = [1,8, 1,9, 6209,120;
    2,7, 2,8, 6453,-240];
end

transX_incremental(1,1) = 0;
transY_incremental(1,1) = 0;
transX(1,1) = 0;
transY(1,1) = 0;



iStop = eNumber + 1;
for iR = 1:size(FOV,1)
    for iC = 1:size(FOV,2)-1
        pause(1);
        close all;
        fName1 = [f1,STOP{iStop},f2,FOV{iR,iC}]
        fName2 = [f1,STOP{iStop},f2,FOV{iR,iC+1}]
        I = imread([path_DIC,'\',fName1,'.tif']);
        J = imread([path_DIC,'\',fName2,'.tif']);
        
        % For debug
        if DRAWFIGURE > 1
            figure;imshowpair(I,J,'montage');
        end
        
        % initially crop a small region to detect
        xi = floor(size(I,2) * 0.01);
        if iR == 1
            yi = floor(size(I,1) * 0.35);
        else
            
            yi = floor(size(I,1) * 0.01);
        end
        
        xf = xi + floor(size(I,2)*0.1);
        yf = yi + floor(size(I,1)*0.1);
        
        Jprime = J(yi:yf,xi:xf);
        if DRAWFIGURE > 0
            figure;imshowpair(I(1:20:end,1:20:end),J(1:20:end,1:20:end),'montage');
        end
        if DRAWFIGURE > 1
            figure;imshowpair(I,Jprime,'montage');
        end
        
        Ipool = my_pool(I,poolSize,1);
        Jpool = my_pool(Jprime,poolSize,1);
        if DRAWFIGURE > 0
            figure;imshowpair(Ipool,Jpool,'montage');
        end
        
        c = normxcorr2(Jpool, Ipool);
        
        if DRAWFIGURE > 1
            figure;surf(c);shading flat; set(gca,'ydir','reverse');view(0,90);
        end
        [ypeak, xpeak] = find(c == max(c(:)));
        xPoolOffSet = xpeak-size(Jpool,2);
        yPoolOffSet = ypeak-size(Jpool,1);
        
        if DRAWFIGURE > 0
            figure;imagesc(Ipool);
            imrect(gca, [xPoolOffSet+1, yPoolOffSet+1, size(Jpool,2), size(Jpool,1)]);
        end
        
        oxi = size(J,2) - xPoolOffSet*poolSize; % overlap_x_initial_size
        oyi = size(J,1) - yPoolOffSet*poolSize; % overlap_y_initial_size
        
        try
            Joverlap = J(1+delta:oyi-delta,1+delta:oxi-delta);
            
            c = normxcorr2(Joverlap, I);
            if DRAWFIGURE > 1
                figure;surf(c);shading flat; set(gca,'ydir','reverse');view(0,90);
            end
            [ypeak, xpeak] = find(c == max(c(:)));
            xOffSet = xpeak-size(Joverlap,2) - delta;
            yOffSet = ypeak-size(Joverlap,1) - delta;
            
            if DRAWFIGURE > 0
                figure;imshowpair(I,J,'montage');
                imrect(gca, [xOffSet+1, yOffSet+1, size(Joverlap,2)+2*delta, size(Joverlap,1)+2*delta]);
            end
            
        catch
            Joverlap = J(1:1000,1:100);
            c = normxcorr2(Joverlap, I);
            
            if DRAWFIGURE > 1
                figure;surf(c);shading flat; set(gca,'ydir','reverse');view(0,90);
            end
            [ypeak, xpeak] = find(c == max(c(:)));
            xOffSet = xpeak-size(Joverlap,2);
            yOffSet = ypeak-size(Joverlap,1);
            
            if DRAWFIGURE > 0
                figure;imshowpair(I,J,'montage');
                imrect(gca, [xOffSet+1, yOffSet+1, size(Joverlap,2)+2*delta, size(Joverlap,1)+2*delta]);
            end
            
        end
        
        transX_incremental(iR,iC+1) = xOffSet
        transY_incremental(iR,iC+1) = yOffSet
        transX(iR,iC+1) = xOffSet + transX(iR,iC)
        transY(iR,iC+1) = yOffSet + transY(iR,iC)
        
        if (iC==1)&&(iR<size(FOV,1))
            pause(1);
            close all;
            fName1 = [f1,STOP{iStop},f2,FOV{iR,iC}]
            fName3 = [f1,STOP{iStop},f2,FOV{iR+1,iC}]
            I = imread([path_DIC,'\',fName1,'.tif']);
            J = imread([path_DIC,'\',fName3,'.tif']);
            
            % For debug
            if DRAWFIGURE > 1
                figure;imshowpair(I,J,'montage');
            end
            
            % initially crop a small region to detect
            xi = floor(size(I,2) * 0.35);
            yi = floor(size(I,1) * 0.01);
            xf = xi + floor(size(I,2)*0.1);
            yf = yi + floor(size(I,1)*0.1);
            
            Jprime = J(yi:yf,xi:xf);
            if DRAWFIGURE > 0
                figure;imshowpair(I(1:20:end,1:20:end),J(1:20:end,1:20:end),'montage');
            end
            if DRAWFIGURE > 1
                figure;imshowpair(I,Jprime,'montage');
            end
            
            Ipool = my_pool(I,poolSize,1);
            Jpool = my_pool(Jprime,poolSize,1);
            if DRAWFIGURE > 0
                figure;imshowpair(Ipool,Jpool,'montage');
            end
            
            c = normxcorr2(Jpool, Ipool);
            if DRAWFIGURE > 1
                figure;surf(c);shading flat; set(gca,'ydir','reverse');view(0,90);
            end
            [ypeak, xpeak] = find(c == max(c(:)));
            xPoolOffSet = xpeak-size(Jpool,2);
            yPoolOffSet = ypeak-size(Jpool,1);
            
            if DRAWFIGURE > 0
                figure;imagesc(Ipool);
                imrect(gca, [xPoolOffSet+1, yPoolOffSet+1, size(Jpool,2), size(Jpool,1)]);
            end
            
            oxi = size(J,2) - xPoolOffSet*poolSize; % overlap_x_initial. 0 = perfect overlap
            oyi = size(J,1) - yPoolOffSet*poolSize;
            
            Joverlap = J(1+delta:oyi-delta,1+delta:oxi-delta);
            c = normxcorr2(Joverlap, I);
            if DRAWFIGURE > 1
                figure;surf(c);shading flat; set(gca,'ydir','reverse');view(0,90);
            end
            [ypeak, xpeak] = find(c == max(c(:)));
            xOffSet = xpeak-size(Joverlap,2) - delta;
            yOffSet = ypeak-size(Joverlap,1) - delta;
            
            if DRAWFIGURE > 0
                figure;imshowpair(I,J,'montage');
                imrect(gca, [xOffSet+1, yOffSet+1, size(Joverlap,2)+2*delta, size(Joverlap,1)+2*delta]);
            end
            
            transX_incremental(iR+1,iC) = xOffSet
            transY_incremental(iR+1,iC) = yOffSet
            transX(iR+1,iC) = xOffSet + transX(iR,iC)
            transY(iR+1,iC) = yOffSet + transY(iR,iC)
            
        end
        
        % if have specialRC to handle
        if exist('specialRC','var')&&(~isempty(specialRC))
            ind = find((iR==specialRC(:,1))&(iC==specialRC(:,2))&(iR==specialRC(:,3))&(iC+1==specialRC(:,4)));
        else
            ind = [];
        end
        if ~isempty(ind)
            transX_incremental(specialRC(ind,3), specialRC(ind,4)) = specialRC(ind,5)
            transY_incremental(specialRC(ind,3), specialRC(ind,4)) = specialRC(ind,6)
            transX(specialRC(ind,3), specialRC(ind,4)) = specialRC(ind,5) + transX(iR,iC)
            transY(specialRC(ind,3), specialRC(ind,4)) = specialRC(ind,6) + transY(iR,iC)
        end
        
    end
end



save(['translations_stop_',num2str(STOP{iStop})],'transX','transY','transX_incremental','transY_incremental');







