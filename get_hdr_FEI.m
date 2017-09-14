% get meta data, and print it into a hdr file (effectively, a txt file)
%
% chenzhe, 2016-10-31

wantStrip = false;    % if you want the 'StripOffsets' and 'StripByteCounts' info, change this to true
if wantStrip
    maxLength = 10e10;
else
    maxLength = 100;
end

directory = uigetdir('','Select a folder that contains all the tif images');

fileList = dir([directory,'\*.tif']);
for ii = 1:size(fileList,1)
    clear str
    fFolder = fileList(ii).folder;  % redundant
    fName = fileList(ii).name;
    metaDataStruct = imfinfo([fFolder,'\',fName]);
    
    % create a new txt file
    fName = strrep(fName,'.tif','.hdr');
    fid = fopen([fFolder,'\',fName],'w');
    
    
    try
        field_DigitalCamera = metaDataStruct.DigitalCamera;
    catch
    end
    
    try
        field_UnknownTags = metaDataStruct.UnknownTags;
    catch
    end
    
    try
        metaDataStruct = rmfield(metaDataStruct,{'DigitalCamera','UnknownTags'});
    catch
    end
    
    if exist('field_DigitalCamera','var')
        fieldNames = fieldnames(field_DigitalCamera);
        for jj = 1:length(fieldNames)
            metaDataStruct = setfield(metaDataStruct,['DigitalCamera',fieldNames{jj}],getfield(field_DigitalCamera,fieldNames{jj}));
        end
    end
    
    if exist('field_UnknownTags','var')
        fieldNames = fieldnames(field_UnknownTags);
        for jj = 1:length(fieldNames)-1
            metaDataStruct = setfield(metaDataStruct,['UnknownTags',fieldNames{jj}],getfield(field_UnknownTags,fieldNames{jj}));
        end
    end
    
    fieldNames = fieldnames(metaDataStruct);
    str = ['[Image File]',char(13),char(10)];
    for jj = 1:length(fieldNames)-2
        str = [str,fieldNames{jj},'='];
        str_a = getfield(metaDataStruct,fieldNames{jj});
        if ~ischar(str_a)
            if length(str_a)<maxLength
                str_a = num2str(str_a);
            else
                str_a = ['length ',num2str(length(str_a)),' vector'];
            end
        end
        str = [str,str_a,char(13),char(10)];
    end
    str = [str,char(13),char(10)];
    fprintf(fid,'%s',str);
    
    try
        fprintf(fid,'%s',field_UnknownTags.Value);
    catch
        str = [];
        for jj = length(fieldNames)-1:length(fieldNames)
            str = [str,fieldNames{jj},'='];
            str_a = getfield(metaDataStruct,fieldNames{jj});
            if ~ischar(str_a)
                if length(str_a)<maxLength
                    str_a = num2str(str_a);
                else
                    str_a = ['length ',num2str(length(str_a)),' vector'];
                end
            end
            str = [str,str_a,char(13),char(10)];
        end
    end
    fprintf(fid,'%s',str);
    fclose(fid);
end
