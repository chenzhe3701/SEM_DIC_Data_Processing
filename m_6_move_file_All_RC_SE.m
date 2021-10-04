% chenzhe, 2017-04-05
% chenzhe, 2017-05-25 make it move between all three folders.
%
% chenzhe, 2017-09-18.
% (1) This code was used to move files among different folders/subfolders.
% (2) There are 3 types of folders: all together, by row-col number, by
% elongation 
% (3) This modification is to use make_FOV_string, so this code can handle
% files that are named in sequence rather than 'rc' format.
% (4) You can modify the appropriate folder/subfolder/filename format, and
% directions of copying.

pathAll = uigetdir('E:\zhec umich Drive\2021-10-01 Mg4Al_m1 insitu SEM-DIC\SEM Data\All','path all');
pathRC = uigetdir('E:\zhec umich Drive\2021-10-01 Mg4Al_m1 insitu SEM-DIC\SEM Data\RC','path to rc');
pathSE = uigetdir('E:\zhec umich Drive\2021-10-01 Mg4Al_m1 insitu SEM-DIC\SEM Data\SE','path to se');
f1 = 'Mg4Al_s';
f2 = '_';
FORMAT = '.mat';    % or '.tif'
copydirection = [2,1];    % Copy direction [from to].  1=all, 2=RC, 3=SE

B = 1;   % 'B' for 'base', to handle if it's 0/1-based index.  But B=1 for 0-based. B=0 for 1-based.  When iR, iC is used with FOV, transX, ... add this B.
row_start = 0;  % starting # of FOV rows
row_end = 3;
col_start = 0;
col_end = 4;    % ending # of FOV cols
e_start = 0;
e_stop = 3;     % elongation #

% file name format: [f1,STOP{#},'_',FOV{#,#}]
% FOV = make_FOV_string(ri, rf, ci, cf, nDigits, sequence)
% sequence = 'rc','snake',or 'raster'
% Usually, all FOVs are analyzed
FOV = make_FOV_string(abs(B-1), row_end, abs(B-1), col_end, 1, 'rc');   


%% 
for iR = row_start:row_end
    for iC = col_start:col_end
        folderRC = [pathRC,'\','r',num2str(iR),'c',num2str(iC)];    % if necessary, change the format of the folder name ------------------
        mkdir(folderRC);    
        for iE = e_start:e_stop
            folderSE = [pathSE,'\','s',num2str(iE)];                 % if necessary, change the format of the folder name ------------------
            mkdir(folderSE);                
            
%             % A simple code if you have more digits to represent elongation leve.  Disable for now.
%             if iE+1<10
%                 STOP = ['00',num2str(iE+1)];
%             else
%                 STOP = ['0',num2str(iE+1)];
%             end
            
            fNameAll = [f1,num2str(iE),'_',FOV{iR+B,iC+B},FORMAT];
            fNameRC = [f1,num2str(iE),'_',FOV{iR+B,iC+B},FORMAT];
            fNameSE = [f1,num2str(iE),'_',FOV{iR+B,iC+B},FORMAT];
            
            copystr{1} = [pathAll,'\',fNameAll];
            copystr{2} = [folderRC,'\',fNameRC];
            copystr{3} = [folderSE,'\',fNameSE];
            % copy from direction 1 to 2
            copyfile(copystr{copydirection(1)},copystr{copydirection(2)},'f');
        end
        disp([iR,iC,iE]);
    end
end
