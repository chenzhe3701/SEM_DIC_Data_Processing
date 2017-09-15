% gei wo de nv shen
% move Marissa's data to a new folder
% the new folder is organized by FOV
% 
% chenzhe, 2017-04-05
% chenzhe, 2017-04-10 update/organize
%
% chenzhe, add note, 2017-09-15
% (1) The use of this code is deprecated. 
% (2) Historically the image names from Marissa's iFast code did not obey a
% common format.
% (3) Initially, images are stored in a parent folder.  Each subfolder
% contains all images of the same elongation.
% (4) This code moves the images into a target parent folder, with
% subfolders, with the image file names corrected.
% (5) You may need to change the folder/file name prefixs (even formats)

folder_source = uigetdir('d:\Marissa_test_20170409','select source image parent folder');
folder_target = uigetdir('d:\Marissa_test_20170409','select target parent folder');

subfolderNamePrefix = '20170430_ts5Al_02_test_e';     % change !!!
fileNamePrefix = '20170430_ts5Al_02_e';     % change !!!
for iR = 0:3
    for iC = 0:13
        for iE = 0:6
            fNameFrom = [fileNamePrefix,num2str(iE),'_','r',num2str(iR),'c',num2str(iC),'.tif'];
            
            % fNameFrom_a handles the file named by ifast.
            if (iC==0)
                fNameFrom_a = [fileNamePrefix,num2str(iE),'_','r',num2str(iR),'c','.tif'];
            else
                fNameFrom_a = [fileNamePrefix,num2str(iE),'_','r',num2str(iR),'c-',num2str(iC),'.tif'];
            end
            
            fNameTo = [fileNamePrefix,num2str(iE),'_','r',num2str(iR),'c',num2str(iC),'.tif'];
            pathFrom = [folder_source,'\',subfolderNamePrefix,num2str(iE)];
            pathTo = [folder_target,'\',subfolderNamePrefix,num2str(iE)];
            mkdir(pathTo);
            try
                copyfile([pathFrom,'\',fNameFrom],[pathTo,'\',fNameTo],'f');
            catch
                copyfile([pathFrom,'\',fNameFrom_a],[pathTo,'\',fNameTo],'f');
            end
        end
        disp([iR,iC]);
    end
end
