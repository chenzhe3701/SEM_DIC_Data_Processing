% zheshigeiwodenvshenxiede. suiranjiandan, danshiwoxiangweitazuoxieshiqing.
%
% combines the code for step 4 and 5 for Marissa's images.
% chenzhe, 20170412 18:56
%
% chenzhe, 2017-05-09, make it specific for these iamge settings
% chenzhe, 2017-05-10, currently I think this is the best solution.
%
% chenzhe, Now just add notes for using code.
%
% chenzhe, 2017-07-17, modify to use my fft_register function.
% It looks less messy now.  The principle is about the same.
%
% chenzhe, 2017-08-06, make this v3, use make_FOV_string().
% currently it is stitch by col.
% If image is too large, cast to 8bit/divide into 4 images.

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
% FOV = make_FOV_string(ri, rf, ci, cf, nDigits, sequence)
% sequence = 'rc','snake',or 'raster'
FOV = make_FOV_string(row_start, row_end, col_start, col_end, 1, 'rc');

clear transX; clear transY; clear transX_incremental; clear transY_incremental;
transX_incremental(1,1) = 0;
transY_incremental(1,1) = 0;
transX(1,1) = 0;
transY(1,1) = 0;
clear specialRC;    % can define special cases
SHIFT_THRESHOLD = 0.75;
%%
for iE = 0:6        % 'e#' in the file name, i.e., stop/pause #  ----------------------------------------------
    folderName = [subFolderNamePrefix_1,num2str(iE)];
    Oly = OVERLAY;
    for iC = col_start:col_end
        for iR = row_start:row_end - 1
            pause(1);
            close all;
            fName1 = [fileNamePrefix_1,num2str(iE),fileNamePrefix_2,'r',num2str(iR),'c',num2str(iC)]; disp(fName1);     % change this accordingly ---------------------------------
            fName2 = [fileNamePrefix_1,num2str(iE),fileNamePrefix_2,'r',num2str(iR+1),'c',num2str(iC)]; disp(fName2);   % change this accordingly -------------------------------
            I = imread([path_DIC,'\',folderName,'\',fName1,'.tif']);
            J = imread([path_DIC,'\',folderName,'\',fName2,'.tif']);
            stitch_direction = 1;   % 1 for up-down, 2 for left-right
            % For debug
            if DRAWFIGURE > 1
                f1 = figure;imshowpair(I,J,'montage');
            end
            
            [yOffSet,xOffSet] = fft_register(I,J,'d',[2500,0, 200,0], [0,2500, 200,0]);         % change this accordingly, be careful with your choice of parameter ---------------------------
            
            transX_incremental(iR+2,iC+1) = xOffSet
            transY_incremental(iR+2,iC+1) = yOffSet
            transX(iR+2,iC+1) = xOffSet  + transX(iR+1,iC+1);
            transY(iR+2,iC+1) = yOffSet  + transY(iR+1,iC+1);
            
            if (iR==row_start)&&(iC<col_end)    % search the right-side picture as well
                pause(1);
                close all;
                fName1 = [fileNamePrefix_1,num2str(iE),fileNamePrefix_2,'r',num2str(iR),'c',num2str(iC)]        % change this accordingly -------------------------------
                fName3 = [fileNamePrefix_1,num2str(iE),fileNamePrefix_2,'r',num2str(iR),'c',num2str(iC+1)]      % change this accordingly -------------------------------
                I = imread([path_DIC,'\',folderName,'\',fName1,'.tif']);
                J = imread([path_DIC,'\',folderName,'\',fName3,'.tif']);
                stitch_direction = 2;   % 1 for up-down, 2 for left-right
                % For debug
                if DRAWFIGURE > 1
                    figure;imshowpair(I,J,'montage');
                end
                
                [yOffSet,xOffSet] = fft_register(I,J,'r',[200,0, 2500,0], [200,0, 100,2500]);       % change this accordingly, be careful with your choice of parameter ---------------------------
                
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
    
    cutEdge = -1;    % cut edge, vs average in blending
    xi = -min(transX(:));     % overall shift when upper_left is (0,0)
    yi = -min(transY(:));
    xm = max(transX(:));    % maximum shift
    ym = max(transY(:));
    img = zeros(resY+yi+ym,resX+xi+xm) * nan;
    
    for iR = 1:size(FOV,1)
        for iC = 1:size(FOV,2)
            fName = [fileNamePrefix_1,num2str(iE),fileNamePrefix_2,FOV{iR,iC}]; disp(fName);    % change this accordingly ----------------------------------------------------------
            
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
    
    bs = whos('img');
    bs = bs.bytes;
    if bs > (2^32-1)*4
        img = uint8(img/255);
        imwrite(img,[path_target,'\img_full_res_',num2str(iE),'.tif'],'tif');
        imwrite(img(1:2:end,1:2:end),[path_target,'\img_half_res_',num2str(iE),'.tif'],'tif');
        imwrite(img(1:10:end,1:10:end),[path_target,'\img_10_res_',num2str(iE),'.tif'],'tif');

    else
        img = uint16(img);
        imwrite(img,[path_target,'\img_full_res_',num2str(iE),'.tif'],'tif');
        imwrite(img(1:2:end,1:2:end),[path_target,'\img_half_res_',num2str(iE),'.tif'],'tif');
        imwrite(img(1:10:end,1:10:end),[path_target,'\img_10_res_',num2str(iE),'.tif'],'tif');
        
%         imwrite(img(1:end/2,1:end/2),[path_target,'\fullImage_1_',num2str(iE),'.tif'],'tif');
%         imwrite(img(1:end/2,end/2+1:end),[path_target,'\fullImage_2_',num2str(iE),'.tif'],'tif');
%         imwrite(img(end/2+1:end,1:end/2),[path_target,'\fullImage_3_',num2str(iE),'.tif'],'tif');
%         imwrite(img(end/2+1:end,end/2+1:end),[path_target,'\fullImage_4_',num2str(iE),'.tif'],'tif');
    end
    imgR = img(1:reduction:end,1:reduction:end);
    imwrite(imgR,[path_target,'\reducedResolutionImage_ver_',num2str(iE),'.tif'],'tif');
end
clock



