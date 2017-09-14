% 2015-10-23, try to combine 36 or 18 FOVs
%
% Currently, the size of the result matrix is directly designed based on
% observation.  Later, the code can be expanded to automatically decide the
% size.  But it is 'later' ... Ti7Al#B6: fov6: x0 = 17176 was wrong, -> 17254, y=0->-10 ->-8
%
% gei nv shen
% 2017-04-03, to Marissa:
% (x=0, y=0) at upper left corner. ...
%
% 2017-04-26, 2017-05-08, 2017-05-09 update. Maybe need to format better, use r#c#. Vic2D-6 might
% have '0' instead of 'NaN'.  Need to look at it's effect.
%
% 2017-05-26, this method is still useful.

[f,p] = uigetfile('D:\Marissa_test_20170430_renamed_cropped_sequential_rotated\translations_searched_vertical_stop_0.mat','select translation');
load([p,f]);    % load translation data
REFSTOP = 1;    % stop number for reference images, default=1

dic_step = 5;       % step size in DIC analysis
try
    % x0 = xTransFOV
    transX = transX
catch
    transX = [];    % can manual input
end
try
    % y0 = yTransFOV
    transY = transY
catch
    transY = [];
end
transX = dic_step * round(transX/dic_step);
transY = dic_step * round(transY/dic_step);
transX_incremental = dic_step * round(transX_incremental/dic_step);
transY_incremental = dic_step * round(transY_incremental/dic_step);

% Marrisa you may use this?
% FOV = {'r0c0'};   % define the name of FOV
minRNum = 0; % starting # of FOV rows
maxRNum = 3;
minCNum = 0;
maxCNum = 13;  % ending # of FOV cols
for iR = minRNum:maxRNum
    for iC = minCNum:maxCNum
        FOV{iR+minRNum+1,iC+minCNum+1} = ['r',num2str(iR),'c',num2str(iC)];
    end
end
nR_FOV = maxRNum-minRNum+1;
nC_FOV = maxCNum-minCNum+1;

% f1 = 'T5_#7_fov_'; f2='_stop_'; f3='.mat';
% Marissa you may use this?
f1 = '20170430_ts5Al_02_e'; f2='_';

% STOP = {'001','002','003','004','005','006','007','008','009','010','011','012','013','014','015'};
% Marissa you may try+change this
STOP = {'0','1','2','3','4','5','6'};

directory_s = uigetdir('D:\Marissa_test_20170430\2]_20170430_ts5Al_02 tensile test_non_incremental_DIC','choose a single directory that directly contains all the DIC mat files');
directory_n = uigetdir('D:\','choose the new/destination directory');

% dic_size = 21800;   % total width of all images.

% Marissa you may want to try the following instead ------------------------- NOTE NOTE NOTE
ii = 0;     % basically, this is x(1) in each FOV. For Vic2D-6, this should always be 0.
resX = 6144; 
resY = 4096;
dic_size_x = 70000;
dic_size_y = 14500;
dataStartX = ii + min(transX(:));   dic_size_x = max(dic_size_x, dataStartX + max(transX(:))+resX)
dataStartY = ii + min(transY(:));   dic_size_y = max(dic_size_y, dataStartY + max(transY(:))+resY)
[X,Y] = meshgrid(dataStartX:dic_step:dic_size_x, dataStartY:dic_step:dic_size_y);

% use 'normal' stitch sequence, or define a special sequence
STITCH_SEQUENCE = 'normal';
if strcmpi(STITCH_SEQUENCE,'normal')
    % default sequence of stitching
    [CS,RS] = meshgrid(1:maxCNum-minCNum+1,1:maxRNum-minRNum+1);
    RS = reshape(RS',1,[]);
    CS = reshape(CS',1,[]);
    iFov_start = 1;
    iFov_stop = length(RS);
    fovR_start = RS(iFov_start);
    fovC_start = CS(iFov_start);
elseif strcmpi(STITCH_SEQUENCE,'special')
    % these are for using RS,CS, it overwrites fovR/C_start. Define it.
    RS = [2*ones(1,10), 3*ones(1,14),   4*ones(1,14),   ones(1,14), 2*ones(1,4)];
    CS = [5:14,         14:-1:1,        1:14,           14:-1:1,    1:4];
    iFov_start = 1;
    iFov_stop = 38;
    fovR_start = RS(iFov_start);
    fovC_start = CS(iFov_start);
end

% fovR_start = 1;     % range in the 'FOV' that you want to merge.
% fovR_stop = 4;
% fovC_start = 1;
% fovC_stop = 14;

% use avg_filter, in order to reduce size later
USE_AVG_FILTER = 0;
Filter_Size = 1;    % odd number, smooth with avg filter, then reduce matrix size.  This can smooth DIC data.
reduce_step = 1;    % reduce data size, data point step size.

% 'xyuv' for raw DIC data, 'XYUV' for distortion corrected, 'Lagrange' for
% distortion corrected and strain calculated.
VarLetter = 'xyuv';
CUTEDGE = 1;    % cut edge or not?
OVERLAY = 200;     % how many overlay to give

%%
for iStop = 1:length(STOP)
    exx = zeros(size(X))*nan;
    exy = zeros(size(X))*nan;
    eyy = zeros(size(X))*nan;
    sigma = ones(size(X))*nan;
    u = zeros(size(X))*nan;
    v = zeros(size(X))*nan;
    
    for iFov = iFov_start:iFov_stop
        fovR = RS(iFov);
        fovC = CS(iFov);
        ind_intercept2 = 0;
        Oly = OVERLAY - 150 + 300 * (iStop-1)/(length(STOP)-1);
        while(sum(ind_intercept2(:))==0)
            Oly = Oly + 150; disp(Oly);
            % fName = [f1,FOV{fovR,fovC},f2,STOP{iStop},f3];
            fName = [f1,STOP{iStop},f2,FOV{fovR,fovC}];  disp(fName);
            fData = load([directory_s,'\',fName,'.mat']);            
            clear mask;
            if (iStop~=REFSTOP)
                try
                    mask = logical(fData.mask);
                    fData.sigma(~mask) = -1;
                    disp('have mask');
                catch
                end
            end
            
            switch VarLetter
                case 'xyuv'             % this works on the raw DIC data
                    u_conv = fData.u;
                    v_conv = fData.v;
                    exx_conv = fData.exx;
                    exy_conv = fData.exy;
                    eyy_conv = fData.eyy;
                case 'XYUV'             % this works best for the data after we performed distortion correction
                    u_conv = fData.U;
                    v_conv = fData.V;
                    exx_conv = fData.exx;
                    exy_conv = fData.exy;
                    eyy_conv = fData.eyy;
                    fData.x = fData.X;
                    fData.y = fData.Y;
                case 'Lagrange'         % this is the data after we perform strain calculation,
                    u_conv = fData.U;
                    v_conv = fData.V;
                    exx_conv = fData.exx_Lagrange;
                    exy_conv = fData.exy_Lagrange;
                    eyy_conv = fData.eyy_Lagrange;
                    fData.x = fData.X;
                    fData.y = fData.Y;
            end
            
            u_conv(fData.sigma==-1)=NaN;
            v_conv(fData.sigma==-1)=NaN;
            exx_conv(fData.sigma==-1)=NaN;
            exy_conv(fData.sigma==-1)=NaN;
            eyy_conv(fData.sigma==-1)=NaN;
            
            if USE_AVG_FILTER
                avgFilter = ones(Filter_Size)/Filter_Size/Filter_Size;
                u_conv = conv2(fData.u, avgFilter, 'same');
                v_conv = conv2(fData.v, avgFilter, 'same');
                exx_conv = conv2(fData.exx, avgFilter, 'same');
                exy_conv = conv2(fData.exy, avgFilter, 'same');
                eyy_conv = conv2(fData.eyy, avgFilter, 'same');
            end
            
            if CUTEDGE
                if (fovR>fovR_start)&&(fovC>fovC_start)
                    fData.sigma(1:size(exx_conv,1)-round((transY(fovR,fovC)-transY(fovR-1,fovC)+Oly)/dic_step), :) = -1;
                    fData.sigma(:, 1:size(exx_conv,2)-round((transX(fovR,fovC)-transX(fovR,fovC-1)+Oly)/dic_step)) = -1;
                    u_conv(1:size(exx_conv,1)-round((transY(fovR,fovC)-transY(fovR-1,fovC)+Oly)/dic_step), :) = NaN;
                    u_conv(:, 1:size(exx_conv,2)-round((transX(fovR,fovC)-transX(fovR,fovC-1)+Oly)/dic_step)) = NaN;
                    v_conv(1:size(exx_conv,1)-round((transY(fovR,fovC)-transY(fovR-1,fovC)+Oly)/dic_step), :) = NaN;
                    v_conv(:, 1:size(exx_conv,2)-round((transX(fovR,fovC)-transX(fovR,fovC-1)+Oly)/dic_step)) = NaN;
                elseif (fovC>fovC_start)     % fovR = 1, cut some on the left
                    fData.sigma(:, 1:size(exx_conv,2)-round((transX(fovR,fovC)-transX(fovR,fovC-1)+Oly)/dic_step)) = -1;
                    u_conv(:, 1:size(exx_conv,2)-round((transX(fovR,fovC)-transX(fovR,fovC-1)+Oly)/dic_step)) = NaN;
                    v_conv(:, 1:size(exx_conv,2)-round((transX(fovR,fovC)-transX(fovR,fovC-1)+Oly)/dic_step)) = NaN;
                elseif (fovR>fovR_start)     % fovC = 1, cut some on the top
                    fData.sigma(1:size(exx_conv,1)-round((transY(fovR,fovC)-transY(fovR-1,fovC)+Oly)/dic_step), :) = -1;
                    u_conv(1:size(exx_conv,1)-round((transY(fovR,fovC)-transY(fovR-1,fovC)+Oly)/dic_step), :) = NaN;
                    v_conv(1:size(exx_conv,1)-round((transY(fovR,fovC)-transY(fovR-1,fovC)+Oly)/dic_step), :) = NaN;
                end
            end
            
            % find index of meaningful data points in global and local coord
            ind_patch_local = (fData.sigma~=-1);    % Locally meaningful data (simga~=-1)
            [nR,nC] = size(ind_patch_local);
            
            C = find(X(1,:) == fData.x(1)+transX(fovR,fovC));   % find starting point of patch in whole matrix.  Note: this code allows some flexibility of the position (+-1 points, currently).
            R = find(Y(:,1) == fData.y(1)+transY(fovR,fovC));
            
            ind_patch_global = zeros(size(X));
            ind_patch_global(R:R+nR-1,C:C+nC-1)=ind_patch_local;    % globally meaningful data (sigma ~=-1)
            
            ind_patch_global = logical(ind_patch_global);     % convert to logical
            ind_patch_local = logical(ind_patch_local);
            
            % sigma
            sigmaPatch = ones(size(X))*nan;     % initialize
            sigmaPatch(ind_patch_global) = fData.sigma(ind_patch_local);    % copy to global position
            if (fovR==fovR_start)&&(fovC==fovC_start)
                sigma = sigmaPatch;
            end
            
            ind_intercept = ~isnan(sigma) & ~isnan(sigmaPatch);     % ind_intercept labels the Intercept region of the FOVs.
            ind_intercept2 = ~isnan(sigma) & (sigmaPatch>-1);       % ind_intercept2 labels the Intercept region of the valid data.
        end
        
        ind_strain = ind_patch_global & (~ind_intercept);       % meaningful data, and not overlapped data.
        if(fovR==fovR_start)&&(fovC==fovC_start)      % special case, 1st image
            ind_strain = ind_patch_global;
        end
        ind_strain_local = ind_strain(R:R+nR-1,C:C+nC-1);
        
        sigma = nanmean(cat(3,sigma,sigmaPatch),3);
        clear sigmaPatch;
        
        % we need to create a special intercept for u/v stitching, making it 'intercept' only on 1 edge, instead of 2 edges.
        % since we are not using [fData.sigma],[sigmaPatch] again, we just modify it, to reduce the overlap area
        % modify ind_intercept, to calculate U/V shift. Copy the 'local'  value, modify, then put it back.  Shrink the intercept region.
        % ind_intercept_local = ind_intercept(R:R+nR-1,C:C+nC-1);
        % ind_intercept(R:R+nR-1,C:C+nC-1) = ind_intercept_local;
        
        
        % U
        uPatch = zeros(size(X))*nan;
        uPatch(ind_patch_global) = u_conv(ind_patch_local);
        if(fovR==fovR_start)&&(fovC==fovC_start)
            u = uPatch;
        end
        uShift1 = u(ind_intercept2);
        uShift1 = nanmean(uShift1(:));
        uShift2 = uPatch(ind_intercept2);
        uShift2 = nanmean(uShift2(:));
        uPatch = uPatch + uShift1 - uShift2;
        u = nanmean(cat(3,u,uPatch),3);
        clear uPatch;
        
        % V
        vPatch = zeros(size(X))*nan;
        vPatch(ind_patch_global) = v_conv(ind_patch_local);
        if (fovR==fovR_start)&&(fovC==fovC_start)
            v = vPatch;
        end
        vShift1 = v(ind_intercept2);
        vShift1 = nanmean(vShift1(:));
        vShift2 = vPatch(ind_intercept2);
        vShift2 = nanmean(vShift2(:));
        vPatch = vPatch + vShift1 - vShift2;
        v = nanmean(cat(3,v,vPatch),3);
        clear vPatch;
        
        % exx, exy, eyy, simply fill with the new one
        exx(ind_strain) = exx_conv(ind_strain_local);
        exy(ind_strain) = exy_conv(ind_strain_local);
        eyy(ind_strain) = eyy_conv(ind_strain_local);        
    end
    sigma(isnan(sigma)) = -1;
    if reduce_step ~= 1
        % data can be taken in reduced density
        x=X(1:reduce_step:end,1:reduce_step:end);
        y=Y(1:reduce_step:end,1:reduce_step:end);
        u=u(1:reduce_step:end,1:reduce_step:end);
        v=v(1:reduce_step:end,1:reduce_step:end);
        exx=exx(1:reduce_step:end,1:reduce_step:end);
        exy=exy(1:reduce_step:end,1:reduce_step:end);
        eyy=eyy(1:reduce_step:end,1:reduce_step:end);
        sigma=sigma(1:reduce_step:end,1:reduce_step:end);
    else
        x=X;
        y=Y;
    end
    save([directory_n,'\',f2,STOP{iStop}],'x','y','u','v','exx','exy','eyy','sigma');
    disp([directory_n,'\',f2,STOP{iStop},' done']);    
end

