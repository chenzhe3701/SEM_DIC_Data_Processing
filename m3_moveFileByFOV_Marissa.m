% gei wo de nv shen
% move Marissa's data to a new folder
% the new folder is organized by FOV
% 
% chenzhe, 2017-04-05
% chenzhe, 2017-04-10 update/organize

folder1 = uigetdir('d:\Marissa_test_20170331','select source image parent folder');
folder2 = uigetdir('d:\Marissa_test_20170331\byFOV','select target parent folder');

pathPrefix = '20170409_ts5Al_01_e';
namePrefix = '20170409_ts5Al_01_e';

for iR = 0:3
    for iC = 0:13
        % make a folder for that FOV
        pathTo = [folder2,'\','r',num2str(iR),'c',num2str(iC)];
        mkdir(pathTo);
        
        % copy all images of this FOV into folder3
        % for each stop
        
        for iE = 0:5
            fNameFrom = [namePrefix,num2str(iE),'_','r',num2str(iR),'c',num2str(iC),'.tif'];
            
            % fNameFrom_a handles the file named by ifast.
            if (iC==0)
                fNameFrom_a = [namePrefix,num2str(iE),'_','r',num2str(iR),'c','.tif'];
            else
                fNameFrom_a = [namePrefix,num2str(iE),'_','r',num2str(iR),'c-',num2str(iC),'.tif'];
            end
            
            fNameTo = [namePrefix,num2str(iE),'_','r',num2str(iR),'c',num2str(iC),'.tif'];
            pathFrom = [folder1,'\',pathPrefix,num2str(iE)];
            try
                copyfile([pathFrom,'\',fNameFrom],[pathTo,'\',fNameTo],'f');
            catch
                copyfile([pathFrom,'\',fNameFrom_a],[pathTo,'\',fNameTo],'f');
            end
        end
        [iR,iC]
    end
end
