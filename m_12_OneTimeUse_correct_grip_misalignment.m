% chenzhezju 2017-5-24
% method 1:
% (1) calculate rotated value, (2) interp to intermediate
% position, (3) intermediate position can simply be pasted to rotated position
%
% chenzhe, 2017-06-08, update.  Using transformation, estimate xi,xf,yi,yf

% chenzhe, 2017-09-20, add note
% (1) I think this code is still useful, so I will edit it. In the future,
% who knows if there will be a test with sample misaligned and a figurable
% misalignment angle? hehe...
% (2) Whether corrected the exy of Vic2D output, this code should both work
% (3) But after running this code, the exy should be corrected! Hurray!
% (4) I think this code was commented well.  So maybe I should just add the
% above notes.
% (5) The suggestion is that: do not only select the undeformed file.

[f,p] = uigetfile('D:\Marissa_test_20170430_stitched_DIC\xyuv_dic_NonIncremental','select dic files need to be rotated','multiselect','on');
if ~iscell(f)
    f = cellstr(f);
end

dic_step = 5;
missAngle = 2.74/180*pi();  % from apparent to what it should be
M3 = angle2dcm(missAngle,0,0,'zxz');
R = M3';
M = M3(1:2,1:2);
% According to defination, X = U*A.
% U is the apparent, X is the actual.  So A should be active, i.e., R
tform = maketform('affine',R);

%%
for iF = 1:length(f)
    % Method (1): try to rotate Rotated back to Original
    % First time, decide the proper range based on 2nd file.
    if 1==iF
        try
            load([p,f{2}],'x','y','u','v','exx','exy','eyy','sigma','exy_corrected');
        catch
            load([p,f{1}],'x','y','u','v','exx','exy','eyy','sigma','exy_corrected');
        end
        rect = [x(1),x(end,1),x(1,end),x(end,end);
            y(1),y(end,1),y(1,end),y(end,end)];
        [rx,ry] = tformfwd(tform, rect(1,:), rect(2,:));
        xi = mod(x(1),dic_step) + floor(min(rx)/dic_step)*dic_step; % mod, in case x(1)~=0.
        xf = mod(x(1),dic_step) + ceil(max(rx)/dic_step)*dic_step;
        yi = mod(y(1),dic_step) + floor(min(ry)/dic_step)*dic_step; % mod, in case x(1)~=0.
        yf = mod(y(1),dic_step) + ceil(max(ry)/dic_step)*dic_step;
        
        [xR,yR] = meshgrid(xi:dic_step:xf,yi:dic_step:yf);
        
        % convert rotated [xR,yR] into intermediate [xInt,yInt], which is position
        % rotated back from rotated-frame to reference-frame (initial)
        xInt = zeros(size(xR))*nan;
        yInt = zeros(size(xR))*nan;
        for ii = 1:size(xR,1)
            for jj = 1:size(xR,2)
                t = [xR(ii,jj), yR(ii,jj)] * M;
                xInt(ii,jj) = t(1);
                yInt(ii,jj) = t(2);
            end
        end
        
        % get a valid sigmaRI/sigmaR, to reduce the size of (xR, yR).
        
        sigmaInt = interp2(x,y,sigma,xInt,yInt,'nearest');
        ri = find(sum((sigmaInt~=-1)&(~isnan(sigmaInt)),2),1,'first');
        rf = find(sum((sigmaInt~=-1)&(~isnan(sigmaInt)),2),1,'last');
        ci = find(sum((sigmaInt~=-1)&(~isnan(sigmaInt)),1),1,'first');
        cf = find(sum((sigmaInt~=-1)&(~isnan(sigmaInt)),1),1,'last');
        
        xInt = xInt(ri:rf,ci:cf);
        yInt = yInt(ri:rf,ci:cf);
        xR = xR(ri:rf,ci:cf);
        yR = yR(ri:rf,ci:cf);
    end
    load([p,f{iF}],'x','y','u','v','exx','exy','eyy','sigma','exy_corrected');
    sigmaR = interp2(x,y,sigma,xInt,yInt,'nearest');
    
    % calculate the rotated values. i.e, do not change coordinates (x,y), but
    % calculate what the value will be after rotation.
    % The for loop looks ugly, but actually might be fater than cellfun() and
    % arrayfun().
    % RI can stand for: value (rotated), position still in reference frame (initial)
    uRI = zeros(size(u))*nan;
    vRI = zeros(size(u))*nan;
    exxRI = zeros(size(u))*nan;
    exyRI = zeros(size(u))*nan;
    eyyRI = zeros(size(u))*nan;
    
    if exist('exy_corrected','var')&&(1==exy_corrected)
        disp('exy already corrected');
        exy_corrected = 1;
    else
        disp('exy being corrected here');
        exy = -exy;
        exy_corrected = 1;
    end
    
    % value rotated, dimension same as x,y,u,v,eij
    for ii = 1:size(u,1)
        for jj = 1:size(u,2)
            t = M*[u(ii,jj); v(ii,jj)];
            uRI(ii,jj) = t(1);
            vRI(ii,jj) = t(2);
            
            t = M*[exx(ii,jj) exy(ii,jj); exy(ii,jj) eyy(ii,jj)]*M';
            exxRI(ii,jj) = t(1);
            exyRI(ii,jj) = t(2);
            eyyRI(ii,jj) = t(4);
        end
    end
    
    % interpolate the values (xxxRI) at positions (xI, yI).  Because (xxxRI) is
    % already rotated, so the interpolated value will be the correct values in
    % the rotated frame.
    
    uR = interp2(x,y,uRI,xInt,yInt,'nearest');
    vR = interp2(x,y,vRI,xInt,yInt,'nearest');
    exxR = interp2(x,y,exxRI,xInt,yInt,'nearest');
    exyR = interp2(x,y,exyRI,xInt,yInt,'nearest');
    eyyR = interp2(x,y,eyyRI,xInt,yInt,'nearest');
    
    %     save([p,f{iF}],'x','y','u','v','exx','exy','eyy','sigma',...
    %         'xI','yI','uRI','vRI','exxRI','exyRI','eyyRI',...
    %         'xR','yR','uR','vR','exxR','exyR','eyyR','sigmaR','-append');
    save([p,f{iF}],'exy','exy_corrected','uRI','vRI','exxRI','exyRI','eyyRI',...
        'xR','yR','uR','vR','exxR','exyR','eyyR','sigmaR','-append');
    clear exy_corrected;
end
