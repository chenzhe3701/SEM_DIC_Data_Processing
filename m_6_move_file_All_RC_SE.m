% move Marissa's data
% chenzhe, 2017-04-05
% chenzhe, 2017-05-25 make it move between all three folders.

pathAll = uigetdir('D:\','path all');
pathRC = uigetdir('D:\Marissa_test_20170430_renamed_cropped_byFOV_rotated\','path to rc');
pathSE = uigetdir('D:\Marissa_test_20170430_renamed_cropped_sequential_rotated\','path to se');
f1 = '20170430_ts5Al_02_e';
f2 = '_';
FORMAT = '.mat';    % or '.tif'
cpdir = [2 1];    % Copy direction [from to].  1=all, 2=RC, 3=SE

for iR = 0:3
    for iC = 0:13
        mkdir([pathRC,'\','r',num2str(iR),'c',num2str(iC)]);
     
%         if (mod(iR,2)==1)
%             iRC = iR*(5+1) + (5-iC)+1;
%         else
%             iRC = iR*(5+1) + iC+1;
%         end
%         if iRC < 10
%             FOV = ['00',num2str(iRC)];
%         else
%             FOV = ['0',num2str(iRC)];
%         end

        for iE = 0:6
           mkdir([pathSE,'\','e',num2str(iE)]);

%             if iE+1<10
%                 STOP = ['00',num2str(iE+1)];
%             else
%                 STOP = ['0',num2str(iE+1)];
%             end
            
            fNameAll = [f1,num2str(iE),'_','r',num2str(iR),'c',num2str(iC),FORMAT];
            fNameRC = [f1,num2str(iE),'_','r',num2str(iR),'c',num2str(iC),FORMAT];
            fNameSE = [f1,num2str(iE),'_','r',num2str(iR),'c',num2str(iC),FORMAT];
            
            cpstr{1} = [pathAll,'\',fNameAll];
            cpstr{2} = [pathRC,'\','r',num2str(iR),'c',num2str(iC),'\',fNameRC];
            cpstr{3} = [pathSE,'\','e',num2str(iE),'\',fNameSE];
            % just change this depending on your copy direction.
            copyfile(cpstr{cpdir(1)},cpstr{cpdir(2)},'f');
        end
        disp([iR,iC,iE]);
    end
end
