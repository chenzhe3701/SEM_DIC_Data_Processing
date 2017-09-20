% columnIndex = find_variable_column_from_CSV_grain_file(EBSDfilePath2, EBSDfileName2, varList)
% EBSDfile is csv file.  Currently at most 30 variables in EBSDfile
% varList = 1 x n cell.  Each cell is the name of one 'variable'. E.g.,
% {'grain-DI','phi1-r'}
% Zhe Chen 20150804 revised.

function columnIndex = find_variable_column_from_CSV_grain_file(EBSDfilePath2, EBSDfileName2, varList)

nVariable = size(varList,2);

fid = fopen([EBSDfilePath2, EBSDfileName2],'r');
c=textscan(fid,'%s',30,'delimiter',',');
columnNames=c{1,1};
fclose(fid);
columnIndex = zeros(nVariable,1);

for iVariable=1:nVariable
    columnIndex(iVariable) = find(strcmpi(columnNames,varList{iVariable}));
end
display('found variable column index from CSV grain file');
display(datestr(now));
end