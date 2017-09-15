% chenzhe 2017-05-27
% modify to include all subfolders
% modify to confirm image size
%
% chenzhe, add notes, 2017-09-15
% Select a parent folder to crop the images in all the subfolders.

start_point = [0 0];
end_point = [6144 4096];
a = inputdlg({'start_x,y and size_x,y'},'Input',1,{num2str([start_point, end_point])});
a = str2num(a{:});
start_point = [a(1),a(2)];
end_point = [a(3),a(4)];

p = uigetdir('F:\','Choose the image (parent) folder');
p = genpath(p);
[a,b] = regexp(p,'[^;]*');
for ii = 1:length(a)
    folder_name{ii} =  [p(a(ii):b(ii)),filesep];
end

for ii = 1:length(folder_name)
    img_names = dir([folder_name{ii} '*.tif']);
    img_names = struct2cell(img_names);
    img_names = img_names(1,:).';
    
    if ~isempty(img_names)
        save_folder_name = [folder_name{ii},'cropped\'];
        mkdir(save_folder_name);
    end        
    
    for jj=1:size(img_names,1)
        img = imread([folder_name{ii}, img_names{jj}],'tif');
        img = imcrop(img,[start_point, end_point]);
        imwrite(img, [save_folder_name, img_names{jj}],'tif')
    end    
end