% [phi1,phi,phi2,x,y,IQ,CI,Phase,Intensity,Fit] = read_ang(fName)
% read the ang file into matrices.
% 
% chenzhe, 2017-07-13.

function [phi1,phi,phi2,x,y,IQ,CI,Phase,Intensity,Fit] = read_ang(fName)
fid = fopen(fName);
C=textscan(fid,'%s','delimiter','\n');
C=C{1,1};
fclose(fid);

% find how many header lines this file has (If line start with #, it is a header)
nRow = size(C,1);
iHeader = 1;
while strcmp('#',C{iHeader}(1))
    iHeader = iHeader+1;
end
nHeader = iHeader - 1;

% initialize
phi1 = zeros(nRow-nHeader,1);
phi = zeros(nRow-nHeader,1);
phi2 = zeros(nRow-nHeader,1);
x = zeros(nRow-nHeader,1);
y = zeros(nRow-nHeader,1);
IQ = zeros(nRow-nHeader,1);
CI = zeros(nRow-nHeader,1);
Phase = zeros(nRow-nHeader,1);
Intensity = zeros(nRow-nHeader,1);
Fit = zeros(nRow-nHeader,1);

% read each row
for ii = (nHeader+1):nRow
    c = textscan(C{ii},'%f');
    c = c{1};
    phi1(ii-nHeader) = c(1);
    phi(ii-nHeader) = c(2);
    phi2(ii-nHeader) = c(3);
    x(ii-nHeader) = c(4);
    y(ii-nHeader) = c(5);
    IQ(ii-nHeader) = c(6);
    CI(ii-nHeader) = c(7);
    Phase(ii-nHeader) = c(8);
    Intensity(ii-nHeader) = c(9);
    Fit(ii-nHeader) = c(10);
end

nC = length(unique(x(:)));
nR = length(unique(y(:)));
phi1 = transpose(reshape(phi1,nC,nR));
phi = transpose(reshape(phi,nC,nR));
phi2 = transpose(reshape(phi2,nC,nR));
x = transpose(reshape(x,nC,nR));
y = transpose(reshape(y,nC,nR));
IQ = transpose(reshape(IQ,nC,nR));
CI = transpose(reshape(CI,nC,nR));
Phase = transpose(reshape(Phase,nC,nR));
Intensity = transpose(reshape(Intensity,nC,nR));
Fit = transpose(reshape(Fit,nC,nR));
end