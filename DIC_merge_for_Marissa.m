% ????? ????? ??
%
% 2017-04-03, to Marissa:
% (x=0, y=0) at upper left corner.
% x0: shift (in pixel) of each image in x position wrt 'r0c0' 
% y0: shift (in pixel) of each image in x position wrt 'r0c0'
% The raw DIC file name is like this: fName=[f1,FOV{fovR,fovC},f2,STOP{iStop},f3].  Change this to load the corresponding DIC data  
% 
% 2015-10-23, try to combine 36 or 18 FOVs
%
% Currently, the size of the result matrix is directly designed based on
% observation.  Later, the code can be expanded to automatically decide the
% size.  But it is 'later' ...

% Ti7Al#B6: fov6: x0 = 17176 was wrong, -> 17254, y=0->-10 ->-8

try 
    x0 = xTransFOV
catch
    x0 = [0	3490	6904	10278	13784	17254
        62	3528	6884	10328	13862	17286
        94	3522	6900	10362	13870	17298
        130	3576	6928	10376	13914	17244
        158	3528	6974	10372	13808	17342
        190	3568	7048	10466	13870	17372];
end
try
    y0 = yTransFOV
catch
    y0 = [0	-4	-8	-8	-8	-8
        3414	3410	3406	3406	3406	3404
        6884	6880	6878	6878	6878	6876
        10210	10204	10200	10198	10198	10196
        13820	13816	13812	13810	13810	13806
        17396	17390	17386	17384	17384	17380];
end

FOV = {
    '001'   '002'   '003'   '004'   '005'   '006'
    '012'   '011'   '010'   '009'   '008'   '007'
    '013'	'014'	'015'	'016'	'017'	'018'
    '024'	'023'	'022'	'021'	'020'	'019'
    '025'	'026'	'027'	'028'	'029'	'030'
    '036'	'035'	'034'	'033'	'032'	'031'};
% Marrisa you may use this? 
% FOV = {'r0c0', 'r0c1', 'r0c2', 'r0c3', 'r0c4', 'r0c5', 'r0c6', 'r0c7', 'r0c8', 'r0c9', 'r0c10', 'r0c11','r0c12', 'r0c13';
% 'r1c0', 'r1c1', 'r1c2', 'r1c3', 'r1c4', 'r1c5', 'r1c6', 'r1c7', 'r1c8', 'r1c9', 'r1c10', 'r1c11','r1c12', 'r1c13';
% 'r2c0', 'r2c1', 'r2c2', 'r2c3', 'r2c4', 'r2c5', 'r2c6', 'r2c7', 'r2c8', 'r2c9', 'r2c10', 'r2c11','r2c12', 'r2c13';
% 'r3c0', 'r3c1', 'r3c2', 'r3c3', 'r3c4', 'r3c5', 'r3c6', 'r3c7', 'r3c8', 'r3c9', 'r3c10', 'r3c11','r3c12', 'r3c13';
% };


f1 = 'T5_#7_fov_'; f2='_stop_'; f3='.mat';
% Marissa you may use this? 
% f1 = '20170331_ts1Al_02_e'; f2='_';

STOP = {'001','002','003','004','005','006','007'};
% STOP = {'001','002','003','004','005','006','007','008','009','010','011','012','013','014','015'};
% Marissa you may try+change this
% STOP = {'0','1','2','3','4','5','6'};

directory_s = uigetdir('','choose the directory of the DIC mat files');
directory_n = uigetdir('','choose the new/destination directory');

dic_step = 3;       % step size in DIC analysis
dic_size = 21800;   % total width of all images.
% Marissa, this maybe ~12000 x 70000, may (or may not) eat up your computer memory...
[X,Y] = meshgrid(1:dic_step:dic_size);    % you can design the size of the data matrix

% Marissa you may want to try the following instead ------------------------- NOTE NOTE NOTE  
% dic_size_x = 21800;
% dic_size_y = 21800;
% [X,Y] = meshgrid(1:dic_step:dic_size_x, 1:dic_step:dic_size_y); 

fovR_start = 1;     % range in the 'FOV' that you want to merge.
fovR_stop = 6;
fovC_start = 1;
fovC_stop = 6;

Filter_Size = 1;    % odd number, smooth with avg filter, then reduce matrix size.  This can smooth DIC data.
reduce_step = 1;    % reduce data size, data point step size.
%%%%%%

%%
for iStop = length(STOP):-1:1
    exx = zeros(size(X))*nan;
    exy = zeros(size(X))*nan;
    eyy = zeros(size(X))*nan;
    sigma = ones(size(X))*nan;
    u = zeros(size(X))*nan;
    v = zeros(size(X))*nan;
    
    for fovR = fovR_start:fovR_stop
        for fovC = fovC_start:fovC_stop
            fName = [f1,FOV{fovR,fovC},f2,STOP{iStop},f3];
            % Marrisa try this?
            % fName = [f1,STOP{iStop},f2,FOV{fovR,fovC}];
            fData = load([directory_s,'\',fName]);
  
            % use avg_filter, in order to reduce size later
            avgFilter = ones(Filter_Size)/Filter_Size/Filter_Size;
            u_conv = conv2(fData.u, avgFilter, 'same');
            v_conv = conv2(fData.v, avgFilter, 'same');
            exx_conv = conv2(fData.exx, avgFilter, 'same');
            exy_conv = conv2(fData.exy, avgFilter, 'same');
            eyy_conv = conv2(fData.eyy, avgFilter, 'same');
            
            % find index of meaningful data points in global and local coord
            ind_patch_local = (fData.sigma~=-1);    % choose data within this patch
            [nR,nC] = size(ind_patch_local);
            
            % On 2016-2-12 I can't remember why I did this, but I disable the following one line by chaning 'ones(13)/13^2' to 'ones(1)/1^2'. !!!!!!!! 
            % 2017-04-03 comment: Yes! '13' was to control the filter/window size.  
            % Only when all the filter is within the image area, the conv result=1, so can be floored to 1, so the point in ind_patch_local can be 1.
            % But this can create 'holes' if sigma==-1 is in a point in the middle rather than edge.  
            ind_patch_local = floor(conv2(double(ind_patch_local), ones(1)/1^2, 'same'));     % looks like this is used to 'cut edge', ZheChen noted 2016-2-12

%             ind_patch_local(1:120,:) = 0;
%             ind_patch_local(nR-120:nR,:) = 0;
%             ind_patch_local(:,1:120) = 0;
%             ind_patch_local(:,nC-120:nC) = 0;
            C = find(X(1,:) == fData.x(1)+x0(fovR,fovC));   % find starting point of patch in whole matrix.  Note: this code allows some flexibility of the position (+-1 points, currently).
            if isempty(C)
                C = find(X(1,:)+1 == fData.x(1)+x0(fovR,fovC)); 
            end
            if isempty(C)
                C = find(X(1,:)-1 == fData.x(1)+x0(fovR,fovC)); 
            end
            R = find(Y(:,1) == fData.y(1)+y0(fovR,fovC));
            if isempty(R)
                R = find(Y(:,1)+1 == fData.y(1)+y0(fovR,fovC));
            end
            if isempty(R)
                R = find(Y(:,1)-1 == fData.y(1)+y0(fovR,fovC));
            end
            ind_patch_global = zeros(size(X));
            ind_patch_global(R:R+nR-1,C:C+nC-1)=ind_patch_local;
            
            ind_patch_global = logical(ind_patch_global);     % convert to logical
            ind_patch_local = logical(ind_patch_local);
            
            % sigma
            sigmaPatch = ones(size(X))*nan;     % initialize
            sigmaPatch(ind_patch_global) = fData.sigma(ind_patch_local);    % copy to global position
            if (fovR==fovR_start)&&(fovC==fovC_start)
                sigma = sigmaPatch;
            end
            
            ind_intercept = ~isnan(sigma) & ~isnan(sigmaPatch);
            ind_strain = ind_patch_global & (~ind_intercept);
            if(fovR==fovR_start)&&(fovC==fovC_start)      % special case, 1st image
               ind_strain = ind_patch_global; 
            end
            ind_strain_local = ind_strain(R:R+nR-1,C:C+nC-1);
            
            sigma = nanmean(cat(3,sigma,sigmaPatch),3);
            clear sigmaPatch;
            
            % U
            uPatch = zeros(size(X))*nan;
            uPatch(ind_patch_global) = u_conv(ind_patch_local);
            if(fovR==fovR_start)&&(fovC==fovC_start)
                u = uPatch;
            end
            uShift1 = u(ind_intercept);
            uShift1 = nanmean(uShift1(:));
            uShift2 = uPatch(ind_intercept);
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
            vShift1 = v(ind_intercept);
            vShift1 = nanmean(vShift1(:));
            vShift2 = vPatch(ind_intercept);
            vShift2 = nanmean(vShift2(:));
            vPatch = vPatch + vShift1 - vShift2;
            v = nanmean(cat(3,v,vPatch),3);
            clear vPatch;
            
            % exx, exy, eyy, simply fill with the new one
            exx(ind_strain) = exx_conv(ind_strain_local);
            exy(ind_strain) = exy_conv(ind_strain_local);
            eyy(ind_strain) = eyy_conv(ind_strain_local);
            
%             exxPatch = zeros(size(X))*nan;
%             exxPatch(ind_patch_global) = exx_conv(ind_patch_local);
%             if (fovR==1)&&(fovC==1)
%                 exx = exxPatch;
%             end
%             exx = nanmean(cat(3,exx,exxPatch),3);
%             clear exxPatch;
%             
%             
%             exyPatch = zeros(size(X))*nan;
%             exyPatch(ind_patch_global) = exy_conv(ind_patch_local);
%             if (fovR==1)&&(fovC==1)
%                 exy = exyPatch;
%             end
%             exy = nanmean(cat(3,exy,exyPatch),3);
%             clear exyPatch;
%             
%             eyyPatch = zeros(size(X))*nan;
%             eyyPatch(ind_patch_global) = eyy_conv(ind_patch_local);
%             if (fovR==1)&&(fovC==1)
%                 eyy = eyyPatch;
%             end
%             eyy = nanmean(cat(3,eyy,eyyPatch),3);
%             clear eyyPatch;
        end
    end
    sigma(isnan(sigma)) = -1;
    % data can be taken in reduced density
    x=X(1:reduce_step:end,1:reduce_step:end);
    y=Y(1:reduce_step:end,1:reduce_step:end);
    u=u(1:reduce_step:end,1:reduce_step:end);
    v=v(1:reduce_step:end,1:reduce_step:end);
    exx=exx(1:reduce_step:end,1:reduce_step:end);
    exy=exy(1:reduce_step:end,1:reduce_step:end);
    eyy=eyy(1:reduce_step:end,1:reduce_step:end);
    sigma=sigma(1:reduce_step:end,1:reduce_step:end);
    save([directory_n,'\YouCouldChangeThisName_',f2,STOP{iStop}],'x','y','u','v','exx','exy','eyy','sigma');
    % Marissa you may try this
    % save([directory_n,'\',f1,STOP{iStop}],'x','y','u','v','exx','exy','eyy','sigma');
    
end
