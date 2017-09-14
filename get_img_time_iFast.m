% read time of image taken.  Start with time = 0 second.
% [ZChen note] Works for images taken in the TESCAN, requiring .hdr files
% 
% chenzhe update 2016-10-31
% Also works for FEI.  First, use 'get_hdr_FEI' to generate hdr files.
% Then run this.
%
% chenzhe update 2017-04-10 
% gei nv shen
% now textCell scans using delimitor 'newline'.

folderName = uigetdir('','choose the folder');
hdrNames = dir([folderName '\*.hdr']);
hdrNames = struct2cell(hdrNames);
hdrNames = hdrNames(1,:).';

imgNameCellArray = cell(length(hdrNames),1);
imgTimeArray = zeros(length(hdrNames),1); 

for ii=1:length(hdrNames)
    fileName = hdrNames{ii};
    imgNameCellArray{ii} = strrep(fileName,'-tif.','.');    % deals with TESCAN hdr file format. 
    imgNameCellArray{ii} = strtok(imgNameCellArray{ii},'.');
%     imgNameCellArray{ii} = strtok(imgNameCellArray{ii},'-');
    
    textCell = textread([folderName '\' fileName], '%s','delimiter','\n');
            
    imgDate = textCell{strncmp(textCell,'Date',4)};
    imgDate = imgDate(6:end);
    imgTime = textCell{strncmp(textCell,'Time',4)};
    imgTime = imgTime(6:end);
    imgTimeArray(ii) = datenum([imgDate,' ',imgTime]);
end

startTime = min(imgTimeArray);

t1 = datevec(imgTimeArray);
t0 = datevec(startTime);

for ii = 1:length(imgTimeArray)
    imgTimeArray(ii) = etime(t1(ii,:),t0);
end

timeData = [{'Image Name'}, {'Time'}; imgNameCellArray, num2cell(imgTimeArray)];

xlswrite([folderName '\time_setup.csv'], timeData);
display('Finish writing to csv file')

