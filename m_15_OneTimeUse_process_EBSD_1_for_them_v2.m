% chenzhe, 2017-06-09
%
% script to correct marissa's EBSD data
%
% chenzhe, 2017_07_10
% an average_affine should be used instead of 'projective'
%
% chenzhe, 2017-09-21
% (1) This code aligns the EBSD and DIC data.  After this the euler angle
% has been changed to the sample reference frame.  However, no symmetry was
% considered, so the changed euler angle is just one of the many solutions.
% (2) Historically for Marissa & Purdue. It might still be useful.
% (3) Note currently it doesn't change the euler angles during alignment.
% (4) Historically, the 20170430 test EBSD data provided to me was flipped,
% so this code has a part to rotate EBSD 180 degrees.  Generally, this
% should not be needed.
% (5) There is no type-2 grain data, and there is no correction to the
% grain statistics data in type-2 grain file if the sample is flipped.
% (6) Other features mentioned in (1), (3) and (5) could be added later.

% addChenFunction;

flipCorrection = 1;     % If the EBSDdata provided to you is flipped.  Mainly historic use.
rotationCorrected = 1;  % If you also have your DIC data rotation corrected, and try to align EBSD to the rotation corrected positions
phiSys = [-90, 180, 0]; % Used to transform the reference frame of the Euler angles.  This setting should work in general, but I need to look at your SEM to confirm.

% csv is currently the (1)cropped & (2)grain dialated data.  I'll leave the raw if you don't like my processing.
[phi1,phi,phi2,x,y,IQ,CI,Fit,ID,edge] = import_EBSD_from_grain_file_1;

% flip marissa's EBSD map in each direction, so it's approximately aligned
% with the sample.
% This is historic.  If someone give you EBSD data that's flipped !!!
% This corrects the Euler angle 180 degree, and change the position.
if flipCorrection == 1
    hwnd = warndlg('Note, you have selected that your EBSD was flipped!');
    waitfor(hwnd);
    phi1 = fliplr(flipud(phi1)) + 180;
    phi = fliplr(flipud(phi));
    phi2 = fliplr(flipud(phi2));
    IQ = fliplr(flipud(IQ));
    CI = fliplr(flipud(CI));
    Fit = fliplr(flipud(Fit));
    ID = fliplr(flipud(ID));
    edge = fliplr(flipud(edge));
end

[phi1,phi,phi2]=align_euler_to_sample(phi1,phi,phi2,'none',phiSys(1),phiSys(2),phiSys(3));  % now align euler with sample.

save('EBSD_euler_aligned_with_sample','phi1','phi','phi2','x','y','IQ','CI','Fit','ID','edge');

% These are the control points.  But if sample is flipped, select them after correction of EBSD. -- I think it's easier to select this way.
% This is for the 20170430 sample.
cpEBSD = [1305, 180;
    1315, 955;
    5765, 485;
    5840, 1120];

cpSEM = [7070,1875;
    6945, 11230;
    60850, 5925;
    61650, 14000];

cpSEM_R = [7150, 1580;
    7475, 10880;
    61100, 3030;
    62250, 11060];

%% project it with the undeformed sample reference frame.

% align EBSD position to the un-rotation-corrected map [x,y]
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

if rotationCorrected == 1
    % align EBSD position to the rotation-corrected map [xR,yR]
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
end
