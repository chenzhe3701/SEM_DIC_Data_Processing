% gei nv shen
%
% chenzhe, 2017-04-10
% stich images


path_DIC = uigetdir('D:\Marissa_test_20170410_cropped\20170409_ts5Al_01_e0','select the folder that contains the images')
path_target = uigetdir('','select a target folder');
[f p] = uigetfile('','select the translation file');
load([p,f]);

% file name format: [f1,STOP{#},'_',FOV{#,#}]
reduction_ratio = 10;
f1 = '20170409_ts5Al_01_e5_';        % change this !!!!!!!!!!!!!!!!!

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


%%
method = 1;
reduction = 10;
resX = 6144;
resY = 4096;
xi = -min(transX(:));     % overall shift
yi = -min(transY(:));
xm = max(transX(:));    % maximum shift
ym = max(transY(:));
img = zeros(resY+yi+ym,resX+xi+xm) * nan;

for iR = 1:size(FOV,1)
    for iC = 1:size(FOV,2)
       fName1 = [f1,FOV{iR,iC}]
       I = double(imread([path_DIC,'\',fName1,'.tif']));    % otherwise, nanmean doesn't work, bucause cat makes nan to 0.
       J = img(1+transY(iR,iC)+yi:resY+transY(iR,iC)+yi,1+transX(iR,iC)+xi:resX+transX(iR,iC)+xi);
       img(1+transY(iR,iC)+yi:resY+transY(iR,iC)+yi,1+transX(iR,iC)+xi:resX+transX(iR,iC)+xi) =...
           mean(cat(3,I,J),3,'omitnan');
    end
    
end
img = uint16(img);
imwrite(img,[path_target,'\fullImage.tif']);
imgR = img(1:reduction:end,1:reduction:end);
imwrite(imgR,[path_target,'\reducedResolutionImage.tif']);

%






