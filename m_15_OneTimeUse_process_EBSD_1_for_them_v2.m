% chenzhe, 2017-06-09
%
% script to correct marissa's EBSD data
%
% chenzhe, 2017_07_10
% an average_affine should be used instead of 'projective'

addChenFunction;
% flip marissa's EBSD map in each direction, so it's approximately aligned
% with the sample.

% csv is currently the (1)cropped & (2)grain dialated data.  I'll leave the raw if you don't like my processing. 
[phi1,phi,phi2,x,y,IQ,CI,Fit,ID,edge] = import_EBSD_1_from_csv;

phi1 = fliplr(flipud(phi1)) + 180; % this corrects the Euler angle due to the simple rotation 
phi = fliplr(flipud(phi));
phi2 = fliplr(flipud(phi2));  
IQ = fliplr(flipud(IQ));
CI = fliplr(flipud(CI));
Fit = fliplr(flipud(Fit));
ID = fliplr(flipud(ID));
edge = fliplr(flipud(edge));

[phi1,phi,phi2]=align_euler_to_sample(phi1,phi,phi2,'none',-90,180,0);  % now align euler with sample

save('EBSD_euler_aligned_with_sample','phi1','phi','phi2','x','y','IQ','CI','Fit','ID','edge');

%% project it with the undeformed sample reference frame.

cpEBSD = [1305, 180;
    1315, 955;
    5765, 485;
    5840, 1120];


% align EBSD position to the un-rotation-corrected map [x,y]

cpSEM = [7070,1875;
    6945, 11230;
    60850, 5925;
    61650, 14000];
tform = make_average_affine(cpEBSD,cpSEM);
[xt,yt] = tformfwd(tform,x,y);
% for simplification
stepSize = 5;
xmin = floor(min(xt(:)));   xmin = xmin - rem(xmin,stepSize);
xmax = ceil(max(xt(:)));    xmax = xmax - rem(xmax,stepSize);
ymin = floor(min(yt(:)));   ymin = ymin - rem(ymin,stepSize);
ymax = ceil(max(yt(:)));    ymax = ymax - rem(ymax,stepSize);
[x_EBSD,y_EBSD] = meshgrid(xmin:stepSize:xmax, ymin:stepSize:ymax);
phi1_EBSD = interp_data(x,y,phi1,x_EBSD,y_EBSD,tform,'interp','nearest');
phi_EBSD = interp_data(x,y,phi,x_EBSD,y_EBSD,tform,'interp','nearest');
phi2_EBSD = interp_data(x,y,phi2,x_EBSD,y_EBSD,tform,'interp','nearest');
IQ_EBSD = interp_data(x,y,IQ,x_EBSD,y_EBSD,tform,'interp','nearest');
CI_EBSD = interp_data(x,y,CI,x_EBSD,y_EBSD,tform,'interp','nearest');
Fit_EBSD = interp_data(x,y,Fit,x_EBSD,y_EBSD,tform,'interp','nearest');
ID_EBSD = interp_data(x,y,ID,x_EBSD,y_EBSD,tform,'interp','nearest');
edge_EBSD = interp_data(x,y,edge,x_EBSD,y_EBSD,tform,'interp','nearest');
save('EBSD_position_aligned_to_unrotated_sample','phi1_EBSD','phi_EBSD','phi2_EBSD','x_EBSD','y_EBSD','IQ_EBSD','CI_EBSD','Fit_EBSD','ID_EBSD','edge_EBSD');


% align EBSD position to the rotation-corrected map [xR,yR]
cpSEM_R = [7150, 1580;
    7475, 10880;
    61100, 3030;
    62250, 11060];
tform = make_average_affine(cpEBSD,cpSEM_R);
[xt,yt] = tformfwd(tform,x,y);
% for simplification
stepSize = 5;
xmin = floor(min(xt(:)));   xmin = xmin - rem(xmin,stepSize);
xmax = ceil(max(xt(:)));    xmax = xmax - rem(xmax,stepSize);
ymin = floor(min(yt(:)));   ymin = ymin - rem(ymin,stepSize);
ymax = ceil(max(yt(:)));    ymax = ymax - rem(ymax,stepSize);
[x_EBSD_R,y_EBSD_R] = meshgrid(xmin:stepSize:xmax, ymin:stepSize:ymax);
phi1_EBSD_R = interp_data(x,y,phi1,x_EBSD_R,y_EBSD_R,tform,'interp','nearest');
phi_EBSD_R = interp_data(x,y,phi,x_EBSD_R,y_EBSD_R,tform,'interp','nearest');
phi2_EBSD_R = interp_data(x,y,phi2,x_EBSD_R,y_EBSD_R,tform,'interp','nearest');
IQ_EBSD_R = interp_data(x,y,IQ,x_EBSD_R,y_EBSD_R,tform,'interp','nearest');
CI_EBSD_R = interp_data(x,y,CI,x_EBSD_R,y_EBSD_R,tform,'interp','nearest');
Fit_EBSD_R = interp_data(x,y,Fit,x_EBSD_R,y_EBSD_R,tform,'interp','nearest');
ID_EBSD_R = interp_data(x,y,ID,x_EBSD_R,y_EBSD_R,tform,'interp','nearest');
edge_EBSD_R = interp_data(x,y,edge,x_EBSD_R,y_EBSD_R,tform,'interp','nearest');
save('EBSD_position_aligned_to_rotated_sample','phi1_EBSD_R','phi_EBSD_R','phi2_EBSD_R','x_EBSD_R','y_EBSD_R','IQ_EBSD_R','CI_EBSD_R','Fit_EBSD_R','ID_EBSD_R','edge_EBSD_R');
