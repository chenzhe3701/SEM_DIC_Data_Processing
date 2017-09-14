% zheshigeiwodenvshenxiede. suiranjiandan, danshiwoxiangweitazuoxieshiqing.
%
% combines the code for step 4 and 5 for Marissa's images.
% chenzhe, 20170412 18:56
%
% chenzhe, 2017-05-09, make it specific for these iamge settings
% chenzhe, 2017-05-10, currently I think this is the best solution.
%
% chenzhe, Now just add notes for using code.

% DRAWFIGURE = 0: no figure, 1: some, 2: more
DRAWFIGURE = 0;     

% Path for images and dic data (dic data usually in the same folder as the images)
path_DIC = uigetdir('D:\Marissa_test_20170430_renamed_cropped\','Select parent folder, which contains subfolders, each containing cropped images at an elongation');
path_target = uigetdir('D:\','select a target folder to hold the stitched images and translation data');

% Sub folder name: [subFolderNamePrefix_1,iE], 
% e.g., 20170430_ts5Al_02_test_e0 
subFolderNamePrefix_1 = '20170430_ts5Al_02_test_e';

% File name format: [fileNamePrefix_1,iE,fileNamePrefix_2='_', 'r', iR, 'c', iC]
% e.g., 20170409_ts5Al_01_e4_r0c0
fileNamePrefix_1 = '20170430_ts5Al_02_e';  
fileNamePrefix_2 = '_';

% resolution of images
resX = 6144;
resY = 4096;
% overlay/window size to search and match images
OVERLAY = 200;
% save reduced size image
reduction = 10;

row_start = 0; % starting # of FOV rows
row_end = 3;
col_start = 0;
col_end = 13;  % ending # of FOV cols

% file name format: [f1,STOP{#},'_',FOV{#,#}]
nFov_row_start = 0; % starting # of FOV rows
nFov_row_end = 3;
nFov_col_start = 0;
nFov_col_end = 13;  % ending # of FOV cols
iR = 0;
for ii = nFov_row_start:nFov_row_end
    iR = iR + 1;
    iC = 0;
    for jj = nFov_col_start:nFov_col_end
        iC = iC + 1;
        FOV{iR,iC} = ['r',num2str(ii),'c',num2str(jj)];
    end
end

clear transX; clear transY; clear transX_incremental; clear transY_incremental;
transX_incremental(1,1) = 0;
transY_incremental(1,1) = 0;
transX(1,1) = 0;
transY(1,1) = 0;
clear specialRC;    % can define special cases
%%
for iE = 0:6        % 'e#' in the file name, i.e., stop/pause #  ----------------------------------------------
    folderName = [subFolderNamePrefix_1,num2str(iE)];
    Oly = OVERLAY;    
    for iC = col_start:col_end
        for iR = row_start:row_end - 1        
            pause(1);
            close all;
            fName1 = [fileNamePrefix_1,num2str(iE),fileNamePrefix_2,'r',num2str(iR),'c',num2str(iC)]; disp(fName1);
            fName2 = [fileNamePrefix_1,num2str(iE),fileNamePrefix_2,'r',num2str(iR+1),'c',num2str(iC)]; disp(fName2);
            I = imread([path_DIC,'\',folderName,'\',fName1,'.tif']);
            J = imread([path_DIC,'\',folderName,'\',fName2,'.tif']);
            stitch_direction = 1;   % 1 for up-down, 2 for left-right
            % For debug
            if DRAWFIGURE > 1
                f1 = figure;imshowpair(I,J,'montage');
            end
            
            try
                % initially crop a small region to detect
                yi = size(I,1) - Oly; % xi, yi are zero-based.
                xi = round(size(I,2)/2-Oly/2);
                
                Iprime = imcrop(I,[xi,yi,Oly,Oly]);
                x_neg = size(I,2)-xi;   % Iprime wrt lower right corner of I
                y_neg = size(I,1)-yi;
                
                c = normxcorr2(Iprime, J);
                
                if DRAWFIGURE > 1
                    figure;surf(c);shading flat; set(gca,'ydir','reverse');view(0,90);
                end
                [ypeak, xpeak] = find(c == max(c(:)));
                xOffSet_J = xpeak-size(Iprime,2);
                yOffSet_J = ypeak-size(Iprime,1);
                
                osx_neg = x_neg + xOffSet_J;    % OffSet_Negative, J's upper left corner's offset wrt I's lower-right corner
                osy_neg = y_neg + yOffSet_J;
                
                xOffSet = size(I,2) - osx_neg;  % positive offset, J's upper-left coner's offset wrt I's upper-left corner
                yOffSet = size(I,1) - osy_neg;
                
                if DRAWFIGURE > 0
                    figure;imshowpair(I,J,'montage');
                    imrect(gca, [xi, yi, size(Iprime,2), size(Iprime,1)]);
                    imrect(gca, [size(I,2)+xOffSet_J, yOffSet_J, size(Iprime,2), size(Iprime,1)]);
                end
                
                if (stitch_direction==1) && ( (yOffSet > size(I,1))||(abs(xOffSet)>500) )
                    error();
                end
                
            catch
                % get rectangular region manually
                try close(f1); catch ; end
                f1 = figure;set(f1,'name','select an area to match, prefer on image I (left)');imshowpair(I,J,'montage');
                rect = round(getrect(gcf));
                xi = rect(1); yi = rect(2); OlyX = rect(3); OlyY = rect(4);
                xi = mod(xi-size(I,2),size(I,2));
                OlyX = mod(OlyX-size(I,2),size(I,2));
                yi = mod(yi-size(I,1),size(I,1));
                OlyY = mod(OlyY-size(I,1),size(I,1));
                Iprime = imcrop(I,[xi,yi,OlyX,OlyY]);
                x_neg = size(I,2)-xi;   % Iprime wrt lower right corner of I
                y_neg = size(I,1)-yi;
                
                c = normxcorr2(Iprime, J);
                if DRAWFIGURE > 1
                    figure;surf(c);shading flat; set(gca,'ydir','reverse');view(0,90);
                end
                [ypeak, xpeak] = find(c == max(c(:)));
                xOffSet_J = xpeak-size(Iprime,2);
                yOffSet_J = ypeak-size(Iprime,1);
                
                osx_neg = x_neg + xOffSet_J;    % J's upper left corner's offset wrt I's lower-right corner
                osy_neg = y_neg + yOffSet_J;
                
                xOffSet = size(I,2) - osx_neg;
                yOffSet = size(I,1) - osy_neg;
                
                if DRAWFIGURE > 0
                    figure;imshowpair(I,J,'montage');
                    imrect(gca, [xi, yi, size(Iprime,2), size(Iprime,1)]);
                    imrect(gca, [size(I,2)+xOffSet_J, yOffSet_J, size(Iprime,2), size(Iprime,1)]);
                end

            end
            
            transX_incremental(iR+2,iC+1) = xOffSet 
            transY_incremental(iR+2,iC+1) = yOffSet 
            transX(iR+2,iC+1) = xOffSet  + transX(iR+1,iC+1);
            transY(iR+2,iC+1) = yOffSet  + transY(iR+1,iC+1);
            
            if (iR==row_start)&&(iC<col_end)    % search the right-side picture as well
                pause(1);
                close all;
                fName1 = [fileNamePrefix_1,num2str(iE),fileNamePrefix_2,'r',num2str(iR),'c',num2str(iC)]
                fName3 = [fileNamePrefix_1,num2str(iE),fileNamePrefix_2,'r',num2str(iR),'c',num2str(iC+1)]
                I = imread([path_DIC,'\',folderName,'\',fName1,'.tif']);
                J = imread([path_DIC,'\',folderName,'\',fName3,'.tif']);
                stitch_direction = 2;   % 1 for up-down, 2 for left-right
                % For debug
                if DRAWFIGURE > 1
                    figure;imshowpair(I,J,'montage');
                end
                
                try
                    % initially crop a small region to detect
                    xi = size(I,2)-Oly; % xi, yi are zero-based.
                    
                    if iR == row_start
                        yi = size(I,1)-Oly - 300;   % 300 treat rotated
                    elseif iR == row_end
                        yi = 0+100;
                    else
                        yi = round(size(I,1)/2-Oly/2);
                    end
                                        
                    Iprime = imcrop(I,[xi,yi,Oly,Oly]);
                    x_neg = size(I,2)-xi;   % Iprime wrt lower right corner of I
                    y_neg = size(I,1)-yi;
                    
                    c = normxcorr2(Iprime, J);
                    
                    if DRAWFIGURE > 1
                        figure;surf(c);shading flat; set(gca,'ydir','reverse');view(0,90);
                    end
                    [ypeak, xpeak] = find(c == max(c(:)));
                    xOffSet_J = xpeak-size(Iprime,2);
                    yOffSet_J = ypeak-size(Iprime,1);
                    
                    osx_neg = x_neg + xOffSet_J;    % J's upper left corner's offset wrt I's lower-right corner
                    osy_neg = y_neg + yOffSet_J;
                    
                    xOffSet = size(I,2) - osx_neg;
                    yOffSet = size(I,1) - osy_neg;
                    
                    if DRAWFIGURE > 0
                        figure;imshowpair(I,J,'montage');
                        imrect(gca, [xi, yi, size(Iprime,2), size(Iprime,1)]);
                        imrect(gca, [size(I,2)+xOffSet_J, yOffSet_J, size(Iprime,2), size(Iprime,1)]);
                     end
                    
                    if (stitch_direction==2) && ( (xOffSet > size(I,2))||(abs(yOffSet)>500) )
                        error();
                    end
                    
                catch
                    % get rectangular region manually
                    try close(f1); catch ; end
                    f1 = figure;set(f1,'name','select an area to match, prefer on image I (left)');imshowpair(I,J,'montage');
                    rect = round(getrect(gcf));
                    xi = rect(1); yi = rect(2); OlyX = rect(3); OlyY = rect(4);
                    xi = mod(xi-size(I,2),size(I,2));
                    OlyX = mod(OlyX-size(I,2),size(I,2));
                    yi = mod(yi-size(I,1),size(I,1));
                    OlyY = mod(OlyY-size(I,1),size(I,1));
                    Iprime = imcrop(I,[xi,yi,OlyX,OlyY]);
                    x_neg = size(I,2)-xi;   % Iprime wrt lower right corner of I
                    y_neg = size(I,1)-yi;
                    
                    c = normxcorr2(Iprime, J);
                    if DRAWFIGURE > 1
                        figure;surf(c);shading flat; set(gca,'ydir','reverse');view(0,90);
                    end
                    [ypeak, xpeak] = find(c == max(c(:)));
                    xOffSet_J = xpeak-size(Iprime,2);
                    yOffSet_J = ypeak-size(Iprime,1);
                    
                    osx_neg = x_neg + xOffSet_J;    % J's upper left corner's offset wrt I's lower-right corner
                    osy_neg = y_neg + yOffSet_J;
                    
                    xOffSet = size(I,2) - osx_neg;
                    yOffSet = size(I,1) - osy_neg;
                    
                    if DRAWFIGURE > 0
                        figure;imshowpair(I,J,'montage');
                        imrect(gca, [xi, yi, size(Iprime,2), size(Iprime,1)]);
                        imrect(gca, [size(I,2)+xOffSet_J, yOffSet_J, size(Iprime,2), size(Iprime,1)]);
                    end
                    
                end
                
                transX_incremental(iR+1,iC+2) = xOffSet
                transY_incremental(iR+1,iC+2) = yOffSet
                transX(iR+1,iC+2) = xOffSet + transX(iR+1,iC+1);
                transY(iR+1,iC+2) = yOffSet + transY(iR+1,iC+1);
                
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
    
    save([path_target,'\','translations_searched_vertical_stop_',num2str(iE)],'transX','transY','transX_incremental','transY_incremental');
    
    cutEdge = 1;    % cut edge, vs average in blending
    xi = -min(transX(:));     % overall shift
    yi = -min(transY(:));
    xm = max(transX(:));    % maximum shift
    ym = max(transY(:));
    img = zeros(resY+yi+ym,resX+xi+xm) * nan;
    
    for iR = 1:size(FOV,1)
        for iC = 1:size(FOV,2)
            fName = [fileNamePrefix_1,num2str(iE),fileNamePrefix_2,FOV{iR,iC}]; disp(fName);
            
            I = double(imread([path_DIC,'\',folderName,'\',fName,'.tif']));    % otherwise, nanmean doesn't work, bucause cat makes nan to 0.
            
            if cutEdge
                if (iR>1)&&(iC>1)
                    I(1:size(I,1)-(transY(iR,iC)-transY(iR-1,iC))-Oly, :) = nan;    % cut top
                    I(:, 1:size(I,2)-(transX(iR,iC)-transX(iR,iC-1))-Oly) = nan;    % cut left
                elseif (iC>1)
                    I(:, 1:size(I,2)-(transX(iR,iC)-transX(iR,iC-1))-Oly) = nan;    % iR=1, cut left
                elseif (iR>1)
                    I(1:size(I,1)-(transY(iR,iC)-transY(iR-1,iC))-Oly, :) = nan;    % iC=1, cut top
                end
            end
            
            J = img(1+transY(iR,iC)+yi:resY+transY(iR,iC)+yi,1+transX(iR,iC)+xi:resX+transX(iR,iC)+xi);
            
            img(1+transY(iR,iC)+yi:resY+transY(iR,iC)+yi,1+transX(iR,iC)+xi:resX+transX(iR,iC)+xi) =...
                mean(cat(3,I,J),3,'omitnan');
        end        
    end
    img = uint16(img);
    imwrite(img,[path_target,'\fullImage_ver_',num2str(iE),'.tif']);
    imgR = img(1:reduction:end,1:reduction:end);
    imwrite(imgR,[path_target,'\reducedResolutionImage_ver_',num2str(iE),'.tif']);        
end
clock



