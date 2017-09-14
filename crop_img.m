
start_point = [0 0];
end_point = [4096 4096];
end_point = [3072 2056];
end_point = [6144 4096];
%folder_name='F:\Images';
folder_name=uigetdir('F:\','Choose the image folder');

if ~strcmp(folder_name(end), '\')
    folder_name = [folder_name '\'];
end
save_folder_name = [folder_name,'cropped\'];
mkdir(save_folder_name);

img_names = dir([folder_name '*.tif']);
img_names = struct2cell(img_names);
img_names = img_names(1,:).';

for ii=1:size(img_names,1)
  img = imread([folder_name, img_names{ii}],'tif');
  img = imcrop(img,[start_point, end_point]);
  imwrite(img, [save_folder_name, img_names{ii}],'tif')
end