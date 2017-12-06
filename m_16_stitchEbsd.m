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
% chenzhe, 2017-08-12.  Almost work.  The fine stitch section still needs
% improvement. Basically, the initial guess for unchecked region needs to
% be updated.
%
% chenzhe, 2017-11-27, note: the 'fit' value was normalized after grain
% dialation, to use to generate grain boundary map.



% DRAWFIGURE = 0: no figure, 1: some, 2: more
DRAWFIGURE = 0;

% Path for images and dic data (dic data usually in the same folder as the images)
path_EBSD_batch = uigetdir('D:\UMich Folder\11 ONR project\31 Experiemt data\[EBSD OIM raw]\WE43_T6_C1 5x27_batch UCSB_2017_08_02','Select parent folder of EBSD data, with cleaned fit value');
path_EBSD_batch_raw = uigetdir('D:\UMich Folder\11 ONR project\31 Experiemt data\[EBSD OIM raw]\WE43_T6_C1 5x27_batch UCSB_2017_08_02','Select parent folder of raw EBSD data');
path_target = uigetdir('D:\','select a target folder to hold the stitched images and translation data');

% File name format: [fileNamePrefix_1, 'r', iR, 'c', iC]
% e.g., 20170409_ts5Al_01_e4_r0c0
fileNamePrefix = 'WE43_T6_C1_';
fileNamePrefix1 = 'Mod.ang_';
fileNamePrefix2 = '_Mod';
fileNamePrefix1_raw = '';   % file name for the raw (not with 'fit' cleaned data)
fileNamePrefix2_raw = '';
fmt = '.ang';

% resolution of images
% resX = 4096;
% resY = 4096;
% overlay/window size to search and match images
OVERLAY = 0;
% save reduced size image
reduction = 10;

row_start = 0; % starting # of FOV rows
row_end = 4;
col_start = 0;
col_end = 26;  % ending # of FOV cols

% FOV = make_FOV_string(ri, rf, ci, cf, nDigits, sequence)
% sequence = 'rc','snake',or 'raster'
FOV = make_FOV_string(0, 4, 0, 26, 1, 'rc');

%%  (1) initial search
clear transX; clear transY; clear transX_incremental; clear transY_incremental;
transX_incremental(1,1) = 0;
transY_incremental(1,1) = 0;
transX(1,1) = 0;
transY(1,1) = 0;
clear specialRC;    % can define special cases
SHIFT_THRESHOLD = 0.75;
resX_max = 0;
resY_max = 0;

Oly = OVERLAY;
for iR = row_start:row_end
    for iC = col_start:col_end - 1
        pause(1);
        close all;
        fName1 = [fileNamePrefix,fileNamePrefix1,'r',num2str(iR),'c',num2str(iC),fileNamePrefix2,fmt]; disp(fName1);     % change this accordingly ---------------------------------
        fName2 = [fileNamePrefix,fileNamePrefix1,'r',num2str(iR),'c',num2str(iC+1),fileNamePrefix2,fmt]; disp(fName2);   % change this accordingly -------------------------------
        [phi1,phi,phi2,~,~,IQ1,~,~,~,Fit1] = read_ang([path_EBSD_batch,'\',fName1]);     % [phi1,phi,phi2,x,y,IQ,CI,Phase,Intensity,Fit] = read_ang(fName1);
        [Phi1,Phi,Phi2,~,~,IQ2,~,~,~,Fit2] = read_ang([path_EBSD_batch,'\',fName2]);
        img1 = find_boundary_from_ID_matrix(Fit1);
        img2 = find_boundary_from_ID_matrix(Fit2);
        
%         img1 = cell2mat(arrayfun(@(x,y,z) calculate_misorientation_euler_d([0,0,0],[x,y,z]/pi*180,'hcp'), phi1, phi, phi2,'uniformoutput',false));
%         img2 = cell2mat(arrayfun(@(x,y,z) calculate_misorientation_euler_d([0,0,0],[x,y,z]/pi*180,'hcp'), Phi1, Phi, Phi2,'uniformoutput',false));
        
%         imwrite(uint16(img1/max(img1(:))*65535),[path_target,'\',fName1,'.tif']);
%         imwrite(uint16(img2/max(img2(:))*65535),[path_target,'\',fName2,'.tif']);
        
        resY_max = max(resY_max,size(img1,1));  % update with size of the biggest FOV
        resX_max = max(resX_max,size(img1,2));
        % For debug
        if DRAWFIGURE > 1
            f1 = figure;imshowpair(img1,img2,'montage');
        end
        
        if(1==iR)
            [yOffSet,xOffSet] = normxcorr2A_register(img2,img1,[1,1,1,1], [40,1,1,1],1);
        else
            [yOffSet,xOffSet] = normxcorr2A_register(img2,img1,[1,1,1,1], [1,1,1,1],1);
        end
        
        transX_incremental(iR+1,iC+2) = xOffSet
        transY_incremental(iR+1,iC+2) = yOffSet
        transX(iR+1,iC+2) = xOffSet  + transX(iR+1,iC+1);
        transY(iR+1,iC+2) = yOffSet  + transY(iR+1,iC+1);
        
        if (iC==col_start)&&(iR<row_end)    % search the down-side picture as well
            pause(1);
            close all;
            fName1 = [fileNamePrefix,fileNamePrefix1,'r',num2str(iR),'c',num2str(iC),fileNamePrefix2,fmt]; disp(fName1);        % change this accordingly -------------------------------
            fName2 = [fileNamePrefix,fileNamePrefix1,'r',num2str(iR+1),'c',num2str(iC),fileNamePrefix2,fmt]; disp(fName2);      % change this accordingly -------------------------------
            [phi1,phi,phi2,~,~,IQ1,~,~,~,Fit1] = read_ang([path_EBSD_batch,'\',fName1]);     % [phi1,phi,phi2,x,y,IQ,CI,Phase,Intensity,Fit] = read_ang(fName1);
            [Phi1,Phi,Phi2,~,~,IQ2,~,~,~,Fit2] = read_ang([path_EBSD_batch,'\',fName2]);
            img1 = find_boundary_from_ID_matrix(Fit1);
            img2 = find_boundary_from_ID_matrix(Fit2);
            
%             img1 = cell2mat(arrayfun(@(x,y,z) calculate_misorientation_euler_d([0,0,0],[x,y,z]/pi*180,'hcp'), phi1, phi, phi2,'uniformoutput',false));
%             img2 = cell2mat(arrayfun(@(x,y,z) calculate_misorientation_euler_d([0,0,0],[x,y,z]/pi*180,'hcp'), Phi1, Phi, Phi2,'uniformoutput',false));
            
            resY_max = max(resY_max,size(img1,1));  % update with size of the biggest FOV
            resX_max = max(resX_max,size(img1,2));
            
            % For debug
            if DRAWFIGURE > 1
                figure;imshowpair(img1,img2,'montage');
            end
            
            if iR==0
                [yOffSet,xOffSet] = normxcorr2A_register(img2,img1,[0,0,0,0], [40,0,0,0], 1);
            else
                [yOffSet,xOffSet] = normxcorr2A_register(img2,img1,[1,1,1,1], [1,1,1,1], 1);
            end
            
            transX_incremental(iR+2,iC+1) = xOffSet
            transY_incremental(iR+2,iC+1) = yOffSet
            transX(iR+2,iC+1) = xOffSet + transX(iR+1,iC+1);
            transY(iR+2,iC+1) = yOffSet + transY(iR+1,iC+1);
            
        end
               
    end
end
disp([resX_max,resY_max]);
save([path_target,'\','translations_ebsd'],'transX','transY','transX_incremental','transY_incremental');
load([path_target,'\','translations_ebsd']);
save([path_target,'\','translations_ebsd_m'],'transX','transY','transX_incremental','transY_incremental');
%% (2) manual confirm.  Can run this multiple times, with different tolerance  
load([path_target,'\','translations_ebsd_m']);

[nR,nC] = size(FOV);    % [nR, nC] is 1-based. right now. need to include 'B' in the future. 

tolerance = 1000;   % if don't want to correct, change to 10000000

buf = 100;          % makes blank space at the edges of the map
xi = -min(transX(:));     % overall shift, leave 100 buffer positions for fine adjustment
yi = -min(transY(:));
xm = max(transX(:));    % maximum shift
ym = max(transY(:));
img = zeros(resY_max+yi+ym + 2*buf, resX_max+xi+xm +2*buf);

transX = transX + buf;    % copy old for initial guess
transY = transY + buf;

% specify an [ir1,ic1; ir2,ic2; ir3,ic3;] to modify
targetRC=[-1 -1];

needInput = 0;
reply = [0 0];
for iRC = 2:nR+nC
    for iR = 1:(iRC-1)
        iC = iRC-iR;
        if (iR<=nR)&&(iC<=nC)
            close all;
            fName = [fileNamePrefix,fileNamePrefix1,FOV{iR,iC},fileNamePrefix2]; disp(fName);
            [~,~,~,~,~,~,~,~,~,Fit] = read_ang([path_EBSD_batch,'\',fName,fmt]);
            img1 = find_boundary_from_ID_matrix(Fit);
            img1(1:10,1:10)=2;
            
            [resY,resX] = size(img1);   % used for copy
            
            img2 = img(1+transY(iR,iC)+yi:resY+transY(iR,iC)+yi,1+transX(iR,iC)+xi:resX+transX(iR,iC)+xi);
            
            img_temp = img;
            img_temp(1+transY(iR,iC)+yi:resY+transY(iR,iC)+yi,1+transX(iR,iC)+xi:resX+transX(iR,iC)+xi) =...
                mean(cat(3,img1,img2),3,'omitnan');
            f=myplot(img_temp);
            set(f,'units','normalized','outerposition',[0 0 1 1]);
            disp([iR iC transX_incremental(iR,iC) transY_incremental(iR,iC)]);
             
            if (abs(transX_incremental(iR,iC) - median(transX_incremental(:,iC)))>tolerance)...
                    || (abs(transY_incremental(iR,iC) - median(transY_incremental(iR,2:end)))>tolerance)...
                    || sum(ismember([iR,iC],targetRC,'rows'))...
                needInput = 1;
                reply = input('input adjustment [x,y]: ','s');
            end
            while (~isempty(reply))&&(needInput)
                reply = str2num(reply);
                transX_incremental(iR,iC) = transX_incremental(iR,iC) + reply(1);
                transY_incremental(iR,iC) = transY_incremental(iR,iC) + reply(2);
                disp([iR iC transX_incremental(iR,iC) transY_incremental(iR,iC)]);
                % update transX transY
                for ir = iR:nR
                    for ic = iC:nC
                        % because searched by row, the same row will be affected 
                        if (ir==iR)
                            transX(ir,ic) = transX(ir,ic) + reply(1);
                        end
                        % if the head column, the transY of the others will be affected  
                        if (1 == iC)||(ir==iR)
                            transY(ir,ic) = transY(ir,ic) + reply(2);
                        end
                    end
                end
                close all;
                
                img2 = img(1+transY(iR,iC)+yi:resY+transY(iR,iC)+yi,1+transX(iR,iC)+xi:resX+transX(iR,iC)+xi);
                
                img_temp = img;
                img_temp(1+transY(iR,iC)+yi:resY+transY(iR,iC)+yi,1+transX(iR,iC)+xi:resX+transX(iR,iC)+xi) =...
                    mean(cat(3,img1,img2),3,'omitnan');
                f = myplot(img_temp);
                set(f,'units','normalized','outerposition',[0 0 1 1]);
                reply = input('input adjustment [x,y]: ','s');
            end
            needInput = 0;
            img = img_temp;
        end
    end
end
myplot(img);
save([path_target,'\','translations_ebsd_m'],'transX','transY','transX_incremental','transY_incremental');

%% (2.0) plot manually adjusted to view
load([path_target,'\','translations_ebsd_m']);
[nR,nC] = size(FOV);

buf = 100;          % makes blank space at the edges of the map
xi = -min(transX(:));     % overall shift, leave 100 buffer positions for fine adjustment
yi = -min(transY(:));
xm = max(transX(:));    % maximum shift
ym = max(transY(:));
img = zeros(resY_max+yi+ym + 2*buf, resX_max+xi+xm +2*buf);
transX = transX + buf;    % copy old for initial guess
transY = transY + buf;

for iRC = 2:nR+nC
    for iR = 1:(iRC-1)
        iC = iRC-iR;
        if (iR<=nR)&&(iC<=nC)
            close all;
            fName = [fileNamePrefix,fileNamePrefix1,FOV{iR,iC},fileNamePrefix2]; disp(fName);
            [~,~,~,~,~,~,~,~,~,Fit] = read_ang([path_EBSD_batch,'\',fName,fmt]);
            img1 = find_boundary_from_ID_matrix(Fit);
            img1(1:10,1:10)=2;
            
            [resY,resX] = size(img1);   % used for copy
            
            img2 = img(1+transY(iR,iC)+yi:resY+transY(iR,iC)+yi,1+transX(iR,iC)+xi:resX+transX(iR,iC)+xi);
            
            img(1+transY(iR,iC)+yi:resY+transY(iR,iC)+yi,1+transX(iR,iC)+xi:resX+transX(iR,iC)+xi) =...
                mean(cat(3,img1,img2),3,'omitnan');
         end
    end
end
myplot(img);


%% (2.1) interp the position.  This should be better.  It uses the adjusted position, but also 'smooth' the adjustment by interp2.
load([path_target,'\','translations_ebsd_m']);
[nR,nC] = size(FOV);

[xq,yq] = meshgrid(1:nC,1:nR);
indR = [1:nR; 1:nR; 1:nR]';
indC = [ones(nR,1)*1, ones(nR,1)*2,ones(nR,1)*nC];
ind = sub2ind([nR,nC],indR,indC);
xx = indC;
yy = indR;
tX = transX(ind);
tY = transY(ind);
transX = round(interp2(xx,yy,tX,xq,yq));
transY = round(interp2(xx,yy,tY,xq,yq));

buf = 100;
xi = -min(transX(:));     % overall shift, leave 100 buffer positions for fine adjustment
yi = -min(transY(:));
xm = max(transX(:));    % maximum shift
ym = max(transY(:));
img = zeros(resY_max+yi+ym + 2*buf, resX_max+xi+xm +2*buf);
transX = transX + buf;    % copy old for initial guess
transY = transY + buf;

for iR = 1:nR
   for iC = 1:nC
      if (1==iR)&&(1==iC)
          transX_incremental(iR,iC) = 0;
          transY_incremental(iR,iC) = 0;
      elseif (iR>1)&&(1==iC)
          % look at top neighbor
          transX_incremental(iR,iC) = transX(iR,iC)-transX(iR-1,iC);
          transY_incremental(iR,iC) = transY(iR,iC)-transY(iR-1,iC);
      elseif (iC>1)
          % look at left neighbor
          transX_incremental(iR,iC) = transX(iR,iC)-transX(iR,iC-1);
          transY_incremental(iR,iC) = transY(iR,iC)-transY(iR,iC-1);
      end
   end
end
for iRC = 2:nR+nC
    for iR = 1:(iRC-1)
        iC = iRC-iR;
        if (iR<=nR)&&(iC<=nC)
            close all;
            fName = [fileNamePrefix,fileNamePrefix1,FOV{iR,iC},fileNamePrefix2]; disp(fName);
            [~,~,~,~,~,~,~,~,~,Fit] = read_ang([path_EBSD_batch,'\',fName,fmt]);
            img1 = find_boundary_from_ID_matrix(Fit);
            img1(1:10,1:10)=2;
            
            [resY,resX] = size(img1);   % used for copy
            
            img2 = img(1+transY(iR,iC)+yi:resY+transY(iR,iC)+yi,1+transX(iR,iC)+xi:resX+transX(iR,iC)+xi);
            
            img(1+transY(iR,iC)+yi:resY+transY(iR,iC)+yi,1+transX(iR,iC)+xi:resX+transX(iR,iC)+xi) =...
                mean(cat(3,img1,img2),3,'omitnan');
         end
    end
end
myplot(img);
save([path_target,'\','translations_ebsd_interp'],'transX','transY','transX_incremental','transY_incremental');

%% (2.2) This is the averaged, like using a constant translation.  Basically, this shows that the stage control can introduce some unreliable drift.
load([path_target,'\','translations_ebsd_m']);
transX_incremental(:,1) = round(mean(transX_incremental(:,1)));
transX_incremental(:,2) = round(mean(transX_incremental(:,2)));
transX_incremental(:,3:end) = round(mean(mean(transX_incremental(:,3:end))));
transY_incremental(:,2:end) = round(mean(mean(transY_incremental(:,2:end))));

buf = 100;
xi = -min(transX(:));     % overall shift, leave 100 buffer positions for fine adjustment
yi = -min(transY(:));
xm = max(transX(:));    % maximum shift
ym = max(transY(:));
img = zeros(resY_max+yi+ym + 2*buf, resX_max+xi+xm +2*buf);
transX = transX + buf;    % copy old for initial guess
transY = transY + buf;

% incremental to transXY
for iR = 1:nR
   for iC = 1:nC
      if (1==iR)&&(1==iC)
          transX(iR,iC) = transX(iR,iC);
          transY(iR,iC) = transY(iR,iC);
      elseif (1==iC)&&(iR>1)
          % look at up
          transX(iR,iC) = transX(iR-1,iC) + transX_incremental(iR,iC);
          transY(iR,iC) = transY(iR-1,iC) + transY_incremental(iR,iC);
      elseif (iC>1)
          % look at left
          transX(iR,iC) = transX(iR,iC-1) + transX_incremental(iR,iC);
          transY(iR,iC) = transY(iR,iC-1) + transY_incremental(iR,iC);
      end       
   end    
end
% plot
for iRC = 2:nR+nC
    for iR = 1:(iRC-1)
        iC = iRC-iR;
        if (iR<=nR)&&(iC<=nC)
            close all;
            fName = [fileNamePrefix,fileNamePrefix1,FOV{iR,iC},fileNamePrefix2]; disp(fName);
            [~,~,~,~,~,~,~,~,~,Fit] = read_ang([path_EBSD_batch,'\',fName,fmt]);
            img1 = find_boundary_from_ID_matrix(Fit);
            img1(1:10,1:10)=2;
            
            [resY,resX] = size(img1);   % used for copy
            
            img2 = img(1+transY(iR,iC)+yi:resY+transY(iR,iC)+yi,1+transX(iR,iC)+xi:resX+transX(iR,iC)+xi);
            
            img(1+transY(iR,iC)+yi:resY+transY(iR,iC)+yi,1+transX(iR,iC)+xi:resX+transX(iR,iC)+xi) =...
                mean(cat(3,img1,img2),3,'omitnan');
         end
    end
end
myplot(img);
save([path_target,'\','translations_ebsd_just_avg'],'transX','transY','transX_incremental','transY_incremental');

%% (3) adjust again for proper size
load([path_target,'\','translations_ebsd_interp']);
transX = transX - min(transX(:));
transY = transY - min(transY(:));
transX_incremental(1,1) = 0;
transY_incremental(1,1) = 0;
for iR = 2:nR
   for iC = 1
       transX_incremental(iR,iC) = transX(iR,iC)-transX(iR-1,iC);
       transY_incremental(iR,iC) = transY(iR,iC)-transY(iR-1,iC);
   end
end
for iR = 1:nR
   for iC = 2:nC
       transX_incremental(iR,iC) = transX(iR,iC)-transX(iR,iC-1);
       transY_incremental(iR,iC) = transY(iR,iC)-transY(iR,iC-1);
   end
end
save([path_target,'\','translations_ebsd_3'],'transX','transY','transX_incremental','transY_incremental');

%% (4) stitch IQ map for illustration, and find actual [nR,nC] of the data
cutEdge = -1;    % cut edge, vs average in blending
xi = -min(transX(:));     % overall shift when upper_left is (0,0)
yi = -min(transY(:));
xm = max(transX(:));    % maximum shift
ym = max(transY(:));
img = zeros(resY_max+yi+ym,resX_max+xi+xm) * nan;
for iR = 1:size(FOV,1)
    for iC = 1:size(FOV,2)
        fName = [fileNamePrefix,fileNamePrefix1,FOV{iR,iC},fileNamePrefix2]; disp(fName);    % change this accordingly ----------------------------------------------------------
        
        [~,~,~,~,~,img1,~,~,~,~] = read_ang([path_EBSD_batch,'\',fName,fmt]);   %[phi1,phi,phi2,x,y,IQ,CI,Phase,Intensity,Fit] = read_ang(fName)
        
        [resY,resX] = size(img1);   % used for copy
        
        img2 = img(1+transY(iR,iC)+yi:resY+transY(iR,iC)+yi,1+transX(iR,iC)+xi:resX+transX(iR,iC)+xi);
        
        img(1+transY(iR,iC)+yi:resY+transY(iR,iC)+yi,1+transX(iR,iC)+xi:resX+transX(iR,iC)+xi) =...
            mean(cat(3,img1,img2),3,'omitnan');
    end
end
data_nR = find(sum(isnan(img),2)==size(img,2),1,'first');
data_nC = find(sum(isnan(img),1)==size(img,1),1,'first');
img = mat2gray(img(1:data_nR,1:data_nC));
imwrite(img,[path_target,'\stitched_IQ','.tif']);

%% (5) stitch real data
phi1 = zeros(data_nR,data_nC) * nan;
phi = phi1;
phi2 = phi1;
x = phi1;
y = phi1;
IQ = phi1;
CI = phi1;
Phase = phi1;
Intensity = phi1;
Fit = phi1;


for iR = 1:size(FOV,1)
    for iC = 1:size(FOV,2)
        fName = [fileNamePrefix,fileNamePrefix1_raw,FOV{iR,iC},fileNamePrefix2_raw]; disp(fName);    % change this accordingly ----------------------------------------------------------
        
        [phi1_t,phi_t,phi2_t,x_t,y_t,IQ_t,CI_t,Phase_t,Intensity_t,Fit_t] = read_ang([path_EBSD_batch_raw,'\',fName,fmt]);       %%[phi1,phi,phi2,x,y,IQ,CI,Phase,Intensity,Fit] = read_ang(fName)
        if (iR==1)&&(iC==1)
            ebsdStep = x_t(1,2)-x_t(1,1);
            [x,y] = meshgrid(0:size(x,2)-1,0:size(x,1)-1);
            x = x*ebsdStep;
            y = y*ebsdStep;
        end
        
        [resY,resX] = size(IQ_t);   % used for copy
        
        IQ_2 = IQ(1+transY(iR,iC)+yi:resY+transY(iR,iC)+yi,1+transX(iR,iC)+xi:resX+transX(iR,iC)+xi);        
        IQ(1+transY(iR,iC)+yi:resY+transY(iR,iC)+yi,1+transX(iR,iC)+xi:resX+transX(iR,iC)+xi) =  mean(cat(3,IQ_t,IQ_2),3,'omitnan');
                
        CI_2 = CI(1+transY(iR,iC)+yi:resY+transY(iR,iC)+yi,1+transX(iR,iC)+xi:resX+transX(iR,iC)+xi);        
        CI(1+transY(iR,iC)+yi:resY+transY(iR,iC)+yi,1+transX(iR,iC)+xi:resX+transX(iR,iC)+xi) =  mean(cat(3, CI_t,CI_2),3,'omitnan');
                
        Phase_2 = Phase(1+transY(iR,iC)+yi:resY+transY(iR,iC)+yi,1+transX(iR,iC)+xi:resX+transX(iR,iC)+xi);        
        Phase(1+transY(iR,iC)+yi:resY+transY(iR,iC)+yi,1+transX(iR,iC)+xi:resX+transX(iR,iC)+xi) =  max(cat(3,Phase_t,Phase_2),[],3);
                
        Intensity_2 = Intensity(1+transY(iR,iC)+yi:resY+transY(iR,iC)+yi,1+transX(iR,iC)+xi:resX+transX(iR,iC)+xi);        
        Intensity(1+transY(iR,iC)+yi:resY+transY(iR,iC)+yi,1+transX(iR,iC)+xi:resX+transX(iR,iC)+xi) =  mean(cat(3,Intensity_t,Intensity_2),3,'omitnan');
                
        Fit_2 = Fit(1+transY(iR,iC)+yi:resY+transY(iR,iC)+yi,1+transX(iR,iC)+xi:resX+transX(iR,iC)+xi);        
        Fit(1+transY(iR,iC)+yi:resY+transY(iR,iC)+yi,1+transX(iR,iC)+xi:resX+transX(iR,iC)+xi) =  mean(cat(3,Fit_t,Fit_2),3,'omitnan');
        
        phi1(1+transY(iR,iC)+yi:resY+transY(iR,iC)+yi,1+transX(iR,iC)+xi:resX+transX(iR,iC)+xi) =  phi1_t;
        phi(1+transY(iR,iC)+yi:resY+transY(iR,iC)+yi,1+transX(iR,iC)+xi:resX+transX(iR,iC)+xi) =  phi_t;
        phi2(1+transY(iR,iC)+yi:resY+transY(iR,iC)+yi,1+transX(iR,iC)+xi:resX+transX(iR,iC)+xi) =  phi2_t;
    end
end

% apply a mask to make the eulers of bad data points = 0.  This helps grain dialation.  
imagesc(IQ); title('Draw polygon over good region. Double click to confirm. Then X to quit.');
H = impoly;
H.wait;
disp('mask created');
mask = H.createMask;
imagesc(mask);title('This is the mask of good data');
mask = double(mask);
phi1(mask==0) = 0;
phi(mask==0) = 0;
phi2(mask==0) = 0;
CI(mask==0) = 0;
IQ(mask==0) = 0;

% remove nan values
phi1(isnan(phi1)) = 0;
phi(isnan(phi)) = 0;
phi2(isnan(phi2)) = 0;
IQ(isnan(IQ)) = 0;
CI(isnan(CI)) = 0;
Phase(isnan(Phase)) = 0;
Intensity(isnan(Intensity)) = 0;
Fit(isnan(Fit)) = 0;

data.phi1 = phi1;
data.phi = phi;
data.phi2 = phi2;
data.x = x;
data.y = y;
data.IQ = IQ;
data.CI = CI;
data.Phase = Phase;
data.Intensity = Intensity;
data.Fit = Fit;
write_ang([path_target,'\stitched_EBSD.ang'],'Mg',data);



