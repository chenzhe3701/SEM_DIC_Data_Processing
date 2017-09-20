% [phi1,phi,phi2,x,y,IQ,CI,Fit,ID,edge] = import_EBSD_1_from_csv(pf)
%
% chenzhe, 2017-06-09
% import EBSD data from csv grain file-1.  Output matrices.

function [phi1,phi,phi2,x,y,IQ,CI,Fit,ID,edge] = import_EBSD_1_from_csv(pf)
    if ~exist('pf','var')
        [f, p] = uigetfile('.csv','choose the EBSD file (csv format, from type-1 grain file)');
        pf=[p,f];
    end
    
    EBSDdata1 = csvread([pf],1,0);
    columnIndex1 = find_variable_column_from_CSV_grain_file([],pf,...
        {'grain-ID','phi1-r','phi-r','phi2-r','x-um','y-um','edge','IQ','CI','Fit'});
    x = EBSDdata1(:,columnIndex1(5));
    y = EBSDdata1(:,columnIndex1(6));
    unique_x = unique(x(:));
    ebsdStepSize = unique_x(2) - unique_x(1);
    mResize = (max(x(:)) - min(x(:)))/ebsdStepSize + 1;
    nResize = (max(y(:)) - min(y(:)))/ebsdStepSize + 1;
    
    phi1 = reshape(EBSDdata1(:,columnIndex1(2)),mResize,nResize)';
    phi = reshape(EBSDdata1(:,columnIndex1(3)),mResize,nResize)';
    phi2 = reshape(EBSDdata1(:,columnIndex1(4)),mResize,nResize)';
    % change it to degrees, if necessary
    if max(phi1(:))<7 && max(phi(:))<7 && max(phi2(:))<7
        phi1 = phi1*180/pi();
        phi = phi*180/pi();
        phi2 = phi2* 180/pi();
    end
    x = reshape(EBSDdata1(:,columnIndex1(5)),mResize,nResize)';
    y = reshape(EBSDdata1(:,columnIndex1(6)),mResize,nResize)';
    ID = reshape(EBSDdata1(:,columnIndex1(1)),mResize,nResize)';
    edge = reshape(EBSDdata1(:,columnIndex1(7)),mResize,nResize)';
    
    IQ = reshape(EBSDdata1(:,columnIndex1(8)),mResize,nResize)';
    CI = reshape(EBSDdata1(:,columnIndex1(9)),mResize,nResize)';
    Fit = reshape(EBSDdata1(:,columnIndex1(10)),mResize,nResize)';    
end