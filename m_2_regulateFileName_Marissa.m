% gei wo de nv shen
% move Marissa's data to a new folder
% the new folder is organized by FOV
% 
% chenzhe, 2017-04-05
% chenzhe, 2017-04-10 update/organize

folder_source = uigetdir('d:\Marissa_test_20170409','select source image parent folder');
folder_target = uigetdir('d:\Marissa_test_20170409','select target parent folder');

pathPrefix = '20170409_ts5Al_01_e';     % change !!!
namePrefix = '20170409_ts5Al_01_e';     % change !!!
for iR = 0:3
    for iC = 0:13
        for iE = 0:5
            fNameFrom = [namePrefix,num2str(iE),'_','r',num2str(iR),'c',num2str(iC),'.tif'];
            
            % fNameFrom_a handles the file named by ifast.
            if (iC==0)
                fNameFrom_a = [namePrefix,num2str(iE),'_','r',num2str(iR),'c','.tif'];
            else
                fNameFrom_a = [namePrefix,num2str(iE),'_','r',num2str(iR),'c-',num2str(iC),'.tif'];
            end
            
            fNameTo = [namePrefix,num2str(iE),'_','r',num2str(iR),'c',num2str(iC),'.tif'];
            pathFrom = [folder_source,'\',pathPrefix,num2str(iE)];
            pathTo = [folder_target,'\',pathPrefix,num2str(iE)];
            mkdir(pathTo);
            try
                copyfile([pathFrom,'\',fNameFrom],[pathTo,'\',fNameTo],'f');
            catch
                copyfile([pathFrom,'\',fNameFrom_a],[pathTo,'\',fNameTo],'f');
            end
        end
        [iR,iC]
    end
end
