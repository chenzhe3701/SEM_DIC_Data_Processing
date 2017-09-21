% chenzhe 2017-05-22
%
% get the value of the specific field in the hdr file

function value = get_hdr_field_value(hdrFullFileName,fieldName)
fid = fopen(hdrFullFileName);
textCell = textscan(fid, '%s','delimiter','\n');
fclose(fid);
textCell = textCell{1};
%textCell = textread(hdrFullFileName, '%s','delimiter','\n');

valueCell = textCell{strncmp(textCell,[fieldName,'='],1+length(fieldName))};

value = valueCell(2+length(fieldName):end);
if ~isempty(str2num(value))
    value = str2num(valueCell(2+length(fieldName):end));
end
end