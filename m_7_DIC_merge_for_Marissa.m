% gei nv shen
%
% 2017-04-03, to Marissa:
% (x=0, y=0) at upper left corner.
% x0: shift (in pixel) of each image in x position wrt 'r0c0'
% y0: shift (in pixel) of each image in x position wrt 'r0c0'
% The raw DIC file name is like this: fName=[f1,FOV{fovR,fovC},f2,STOP{iStop},f3].  Change this to load the corresponding DIC data
%
% 2015-10-23, try to combine 36 or 18 FOVs
%
% Currently, the size of the result matrix is directly designed based on
% observation.  Later, the code can be expanded to automatically decide the
% size.  But it is 'later' ...
%
% Ti7Al#B6: fov6: x0 = 17176 was wrong, -> 17254, y=0->-10 ->-8
%
% 2017-04-26, 2017-05-08, 2017-05-09 update. Maybe need to format better, use r#c#. Vic2D-6 might
% have '0' instead of 'NaN'.  Need to look at it's effect.
%
% 2017-05-14, rewrite a lot.
% Now the stitching is based on an least-square fit algorithm, trying ot
% minimize the error of the misfit of the u/v values on the overlap regions
% by provide an additional u/vTrans to the translation data.
%
% However, the stitching still uses the 'cutEdge' method, rather than using
% a similar method in stitching exx/exy/eyy (which just paste data to a new
% regions that has not been filled)
%
% 2017-09-19, chenzhe
% combine the features from v3 and v_NewAndOld, and write this code.
% (1) You can choose to use a global method to calculate the adjustment
% that should be applied to uv displacement values to make the stitched
% displacement fields continuous.
% (2) You can choose to use a local method, that stitch each of the FOV
% one-by-one, and calculate this adjustment everytime.
% (3) If use the local method, you can also define a special sequence to
% stitch each FOVs.

[f,p] = uigetfile('D:\Marissa_test_20170430_renamed_cropped_sequential_rotated\translations_searched_vertical_stop_0.mat','select translation');
load([p,f]);    % load translation data
iE_ref = 0;    % elongation number for reference images, default=0

dic_step = 5;       % step size in DIC analysis

% transX = [];    % can manual input
% transY = [];

transX = dic_step * round(transX/dic_step);
transY = dic_step * round(transY/dic_step);
transX_incremental = dic_step * round(transX_incremental/dic_step);
transY_incremental = dic_step * round(transY_incremental/dic_step);

B = 1;   % 'B' for 'base', to handle if it's 0/1-based index.  But B=1 for 0-based. B=0 for 1-based.  When iR, iC is used with FOV, transX, ... add this B.
row_begin = 0; % starting # of FOV rows
row_end = 3;
col_begin = 0;
col_end = 13;  % ending # of FOV cols
e_begin = 0;
e_end = 6;

FOV = make_FOV_string(abs(B-1), row_end, abs(B-1), col_end, 1, 'rc');

% STOP = {'0','1','2','3','4','5','6'};
for ii = e_begin:e_end
    STOP{ii+B} = num2str(ii);
end

nR_FOV = row_end-row_begin+1;
nC_FOV = col_end-col_begin+1;

stitch_method = 'global';  % (1) global, (2) local, one-by-one
% Can use (1) 'normal' stitch sequence, or (2) define a special sequence 
stitch_sequence = 'normal';
if strcmpi(stitch_sequence,'normal')
    % default sequence of stitching
    [CS,RS] = meshgrid(col_begin:col_end,row_begin:row_end);
    RS = reshape(RS',1,[]);
    CS = reshape(CS',1,[]);
    iFOV_begin = 1;
    iFOV_end = length(RS);
elseif strcmpi(stitch_sequence,'special')&&strcmpi(stitch_method,'local')
    % these are for using RS,CS. define it.
    RS = [1*ones(1,10), 2*ones(1,14),   3*ones(1,14),   0*ones(1,14), 1*ones(1,4)];
    CS = [4:13,         13:-1:0,        0:13,           13:-1:0,    0:3];
    iFOV_begin = 1;
    iFOV_end = 38;
else
    error('global-method with special-sequence not allowed');
end

% File name format: [fileNamePrefix_1,iE,fileNamePrefix_2='_', 'r', iR, 'c', iC]
% e.g., 20170409_ts5Al_01_e4_r0c0
f1 = '20170430_ts5Al_02_e'; f2='_';

directory_s = uigetdir('D:\Marissa_test_20170430\4]_20170430_ts5Al_02 tensile test_non_incremental_DIC_rotated','choose a single directory that directly contains all the DIC mat files');
directory_n = uigetdir('D:\','choose the new/destination directory');

% Basically, this is the smallest x/y position in each FOV.
% For Vic2D-6 this should be 0.
% For Vic2D-2009, with subset=21, step = 5, you could have x=12,17,..., so
% in this case xyi=2.  Then also you should have the same x=12,17,..., for
% all your FOVs.
xyi = 0;

resX = 6144;            % res of each image
resY = 4096;

xi_mesh = xyi + min(transX(:));
xf_mesh = xi_mesh + max(transX(:))+resX + (max(transX(:,1))-xi_mesh);
yi_mesh = xyi + min(transY(:));
yf_mesh = yi_mesh + max(transY(:))+resY + (max(transY(1,:))-yi_mesh);
[X,Y] = meshgrid(xi_mesh:dic_step:xf_mesh+dic_step, yi_mesh:dic_step:yf_mesh+dic_step);

% use avg_filter, in order to reduce size later
USE_AVG_FILTER = 0;
Filter_Size = 1;    % odd number, smooth with avg filter, then reduce matrix size.  This can smooth DIC data.
reduce_step = 1;    % reduce data size, data point step size.

% 'xyuv' for raw DIC data, 'XYUV' for distortion corrected, 'Lagrange' for
% distortion corrected and strain calculated.
VarLetter = 'xyuv';
CUTEDGE = 1;        % cut edge or not?
OVERLAY = 200;      % how many overlay to give

pts_considered = 10;      % how many points in each row/column to determine u/v shift

pts_ignored = 10 + 1;           % ignore some additional (rows or columns) points in the overlap region
pts_ignored_left = 300 + 1;     % for image taken withou external scan, the distortion is larger, so igmore more points


%%
for iE = e_begin:e_end
    exx = zeros(size(X))*nan;
    exy = zeros(size(X))*nan;
    eyy = zeros(size(X))*nan;
    sigma = ones(size(X))*nan;
    u = zeros(size(X))*nan;
    v = zeros(size(X))*nan;
    
    
    if strcmpi(stitch_method,'global')
        % Prepare, to calculate the uvTrans (which is an adjustment) that can be applied to the displacement data.
        N = nR_FOV*nC_FOV;  % total number of elements/unknown translations/FOVs
        
        % E=Sigma(u1+uTrans1-u2-uTrans2)^2
        % 0=u1+uTrans1-u2-uTrans2
        % Y = (K.*A)*X, Y=u2-u1, X=unknown(i.e.,uTrans)
        A = zeros(N^2,N);
        Yu = zeros(N^2,1);
        Yv = zeros(N^2,1);
        K = zeros(N^2,1);
        
        % B1=(u1-u2)^2, B2=2(u1-u2)(uTrans1-uTrans2),B3=(uTrans1-uTrans2)^2
        B1 = zeros(N,N);
        B2 = zeros(1,N);
        B3 = 0;
        D1 = zeros(N,N);
        D2 = zeros(1,N);
        D3 = 0;
        
        for ii = 1:N                % ii = 1st fov #
            for jj = ii+1:N         % jj = 2nd fov #
                rr = (ii-1)*N + jj; % rr = indr in A matrix
                A(rr,ii) = 1;
                A(rr,jj) = -1;
            end
        end
        
        % [iR,iC] are for FOVs.  For each FOV and its neighbor, calculate some coefficients
        for iR = row_begin:row_end
            for iC = col_begin:col_end
                % (0) current fov data
                fName = [f1,STOP{iE+B},f2,FOV{iR+B,iC+B}]; disp(fName);
                fData = load([directory_s,'\',fName,'.mat']);
                clear mask;
                if (iE~=iE_ref)
                    try
                        mask = logical(fData.mask);
                        fData.sigma(~mask) = -1;
                        disp('have mask');
                    catch
                    end
                end
                flag0 = double((fData.sigma~=-1)&(~isnan(fData.sigma)));    % use this to label overlap region by looking at sigma values
                
                % [nR,nC] are for data points. Index of this patch in the global matrix is [Ri:Rf, Ci:Cf]
                [nR,nC] = size(fData.sigma);
                [~,Ci] = min(abs(fData.x(1)+transX(iR+B,iC+B) - X(1,:)));   % find starting point of patch in whole matrix.
                [~,Ri] = min(abs(fData.y(1)+transY(iR+B,iC+B) - Y(:,1)));
                Cf = Ci + nC - 1;
                Rf = Ri + nR - 1;
                
                % (1) search for right neighbor
                if iC < col_end
                    fName_1 = [f1,STOP{iE+B},f2,FOV{iR+B,iC+B+1}]; disp(fName_1);
                    fData_1 = load([directory_s,'\',fName_1,'.mat']);
                    clear mask;
                    if (iE~=iE_ref)
                        try
                            mask = logical(fData_1.mask);
                            fData_1.sigma(~mask) = -1;
                            disp('have mask');
                        catch
                        end
                    end
                    flag_1 = double((fData_1.sigma~=-1)&(~isnan(fData_1.sigma)));    % use this to label overlap region by looking at sigma values
                    
                    % Similarly, this patch in the global matrix is indexed as [Ri_1:Rf_1, Ci_1:Cf_1]
                    [nR_1,nC_1] = size(fData_1.sigma);
                    [~,Ci_1] = min(abs(fData_1.x(1)+transX(iR+B,iC+B+1) - X(1,:)));   % find starting point of patch in whole matrix.
                    [~,Ri_1] = min(abs(fData_1.y(1)+transY(iR+B,iC+B+1) - Y(:,1)));
                    Cf_1 = Ci_1 + nC_1 - 1;
                    Rf_1 = Ri_1 + nR_1 - 1;
                    
                    % for a FOV, index of the region overlapped with its right neighbor
                    ci = max(1, 1+Ci_1-Ci);
                    cf = min(nC, nC + Cf_1-Cf);
                    ri = max(1, 1+Ri_1-Ri);
                    rf = min(nR, nR + Rf_1-Rf);
                    
                    % on the neighbor FOV, index of this overlap region
                    ci_1 = ci + Ci-Ci_1;
                    cf_1 = cf + Ci-Ci_1;
                    ri_1 = ri + Ri-Ri_1;
                    rf_1 = rf + Ri-Ri_1;
                    
                    % Update flag to label the regions (of course, still within overlap regions) on the FOV to be used to calculate uvShift
                    flag = flag0;
                    flag(ri:rf,ci:cf) = flag(ri:rf,ci:cf) + flag_1(ri_1:rf_1,ci_1:cf_1);
                    flag = double(flag==2);       % This is the overlay region
                    iii = find(sum(flag,2),pts_ignored);    % do not consider the 1st or last 10 rows if analyzing with right neighbor (cols if analyzing with down neighbor)
                    iii = iii(end);
                    iij = find(sum(flag,2),pts_ignored,'last');
                    iij = iij(1);
                    for iRR = iii:iij
                        indC = find(flag(iRR,:),pts_considered,'last');
                        if ~isempty(indC)
                            indC = indC(1:min(pts_considered,length(indC)));
                            flag(iRR, indC) = 2;       % for selected rows, label nSamples (if not enough, label all of the) points
                        end
                    end
                    flag = (flag==2); % logical, update flagT with the points just labeled
                    % Similarly, update flag_1
                    flag_1 = false(size(flag_1));
                    flag_1(ri_1:rf_1,ci_1:cf_1) = flag(ri:rf,ci:cf);
                    
                    % copy the measured u,v value of the intereted regions for the two FOVs
                    uu = fData.u(flag);
                    uu_1 = fData_1.u(flag_1);
                    vv = fData.v(flag);
                    vv_1 = fData_1.v(flag_1);
                    
                    % find the ind in the coefficient matrix of the equations to be solved
                    ind = (iR+B-1)*nC_FOV + iC+B;
                    ind_1 = (iR+B-1)*nC_FOV + iC+B+1;
                    rr = (ind-1)*N + ind_1;
                    
                    % calculate terms in the equations
                    K(rr) = K(rr)+length(uu);
                    Yu(rr) = Yu(rr) + sum(uu_1) - sum(uu);
                    Yv(rr) = Yv(rr) + sum(vv_1) - sum(vv);
                    
                    B1(ind,ind) = B1(ind,ind) + length(uu);
                    B1(ind_1,ind_1) = B1(ind_1,ind_1) + length(uu);
                    B1(ind,ind_1) = B1(ind,ind_1) - 2*length(uu);
                    B2(ind) = B2(ind) + 2 * sum(uu - uu_1);
                    B2(ind_1) = B2(ind_1) - 2 * sum(uu - uu_1);
                    B3 = B3 + (uu - uu_1)'*(uu - uu_1);
                    
                    D1(ind,ind) = D1(ind,ind) + length(vv);
                    D1(ind_1,ind_1) = D1(ind_1,ind_1) + length(vv);
                    D1(ind,ind_1) = D1(ind,ind_1) - 2*length(vv);
                    D2(ind) = D2(ind) + 2 * sum(vv - vv_1);
                    D2(ind_1) = D2(ind_1) - 2 * sum(vv - vv_1);
                    D3 = D3 + (vv - vv_1)'*(vv - vv_1);
                end
                
                % (2) search for the bottom neighbor
                if iR < row_end
                    fName_1 = [f1,STOP{iE+B},f2,FOV{iR+B+1,iC+B}]; disp(fName_1);
                    fData_1 = load([directory_s,'\',fName_1,'.mat']);
                    clear mask;
                    if (iE~=iE_ref)
                        try
                            mask = logical(fData_1.mask);
                            fData_1.sigma(~mask) = -1;
                            disp('have mask');
                        catch
                        end
                    end
                    flag_1 = double((fData_1.sigma~=-1)&(~isnan(fData_1.sigma)));    % use this to label overlap region by looking at sigma values
                    
                    [nR_1,nC_1] = size(fData_1.sigma);
                    [~,Ci_1] = min(abs(fData_1.x(1)+transX(iR+B+1,iC+B) - X(1,:)));   % find starting point of patch in whole matrix.
                    [~,Ri_1] = min(abs(fData_1.y(1)+transY(iR+B+1,iC+B) - Y(:,1)));
                    Cf_1 = Ci_1 + nC_1 - 1;
                    Rf_1 = Ri_1 + nR_1 - 1;
                    
                    ci = max(1, 1+Ci_1-Ci);
                    cf = min(nC, nC + Cf_1-Cf);
                    ri = max(1, 1+Ri_1-Ri);
                    rf = min(nR, nR + Rf_1-Rf);
                    
                    ci_1 = ci + Ci-Ci_1;
                    cf_1 = cf + Ci-Ci_1;
                    ri_1 = ri + Ri-Ri_1;
                    rf_1 = rf + Ri-Ri_1;
                    
                    flag = flag0;
                    flag(ri:rf,ci:cf) = flag(ri:rf,ci:cf) + flag_1(ri_1:rf_1,ci_1:cf_1);
                    flag = double(flag==2);   % This is the overlay region
                    iii = find(sum(flag,1),pts_ignored_left);   % throw away left data
                    iii = iii(end);
                    iij = find(sum(flag,1),pts_ignored,'last');
                    iij = iij(1);
                    for iCC = iii:iij
                        indR = find(flag(:,iCC),pts_considered,'last');
                        if ~isempty(indR)
                            indR = indR(1:min(pts_considered,length(indR)));
                            flag(indR,iCC) = 2;
                        end
                    end
                    flag = (flag==2); % logical
                    flag_1 = false(size(flag_1));
                    flag_1(ri_1:rf_1,ci_1:cf_1) = flag(ri:rf,ci:cf);
                    uu = fData.u(flag);
                    uu_1 = fData_1.u(flag_1);
                    vv = fData.v(flag);
                    vv_1 = fData_1.v(flag_1);
                    
                    ind = (iR+B-1)*nC_FOV + iC+B;
                    ind_1 = (iR+B-1+1)*nC_FOV + iC+B;
                    rr = (ind-1)*N + ind_1;
                    
                    K(rr) = K(rr)+length(uu);
                    Yu(rr) = Yu(rr) + sum(uu_1) - sum(uu);
                    Yv(rr) = Yv(rr) + sum(vv_1) - sum(vv);
                    
                    B1(ind,ind) = B1(ind,ind) + length(uu);
                    B1(ind_1,ind_1) = B1(ind_1,ind_1) + length(uu);
                    B1(ind,ind_1) = B1(ind,ind_1) - 2*length(uu);
                    B2(ind) = B2(ind) + 2 * sum(uu - uu_1);
                    B2(ind_1) = B2(ind_1) - 2 * sum(uu - uu_1);
                    B3 = B3 + (uu - uu_1)'*(uu - uu_1);
                    
                    D1(ind,ind) = D1(ind,ind) + length(vv);
                    D1(ind_1,ind_1) = D1(ind_1,ind_1) + length(vv);
                    D1(ind,ind_1) = D1(ind,ind_1) - 2*length(vv);
                    D2(ind) = D2(ind) + 2 * sum(vv - vv_1);
                    D2(ind_1) = D2(ind_1) - 2 * sum(vv - vv_1);
                    D3 = D3 + (vv - vv_1)'*(vv - vv_1);
                end
                disp(['RC: ',num2str(iR),' , ',num2str(iC)]);
            end
        end
        
        % Equation 0.
        ind = (K~=0);
        K = K(ind);
        A = A(ind,:);
        Yu = Yu(ind);
        Yv = Yv(ind);
        KA = K.*A;      % KA is rank deficient, but LS is left inverse.
        
        uTrans = KA\Yu;
        uTrans = transpose(reshape(uTrans-uTrans(1),nC_FOV,nR_FOV));
        vTrans = KA\Yv;
        vTrans = transpose(reshape(vTrans-vTrans(1),nC_FOV,nR_FOV));
        
        % Equation 1. OK. Solution is here, may try it.
        uTrans1 = fminunc(@(x) sum(sum(B1.*(x*x')))+B2*x+B3, reshape(uTrans',[],1));
        uTrans1 = transpose(reshape(uTrans1-uTrans1(1),nC_FOV,nR_FOV));
        vTrans1 = fminunc(@(x) sum(sum(D1.*(x*x')))+D2*x+D3, reshape(uTrans',[],1));
        vTrans1 = transpose(reshape(vTrans1-vTrans1(1),nC_FOV,nR_FOV));
        
        % Equation 2. This doesn't work good
        uTrans2 = inv(KA'*KA)*KA'*Yu;
        uTrans2 = transpose(reshape(uTrans2-uTrans2(1),nC_FOV,nR_FOV));
        vTrans2 = inv(KA'*KA)*KA'*Yv;
        vTrans2 = transpose(reshape(vTrans2-vTrans2(1),nC_FOV,nR_FOV));
        save([directory_n,'\','methodGlobal_uvTrans_e',STOP{iE+B}],'uTrans','vTrans','uTrans1','vTrans1','uTrans2','vTrans2');
    end
    
    
    
    
    % Now, loop through each FOV and stitch data.
    for ii = iFOV_begin:iFOV_end
        iR = RS(ii);
        iC = CS(ii);
        
        % first, find intercept region
        ind_intercept = 0;
        Oly = OVERLAY-150 + 300 * iE/length(STOP);  % ---------- this is actually based on experience ... -------------------- 
        while(sum(ind_intercept(:))==0)
            Oly = Oly + 150; disp(['current overly: ',num2str(Oly)]);
            
            fName = [f1,STOP{iE+B},f2,FOV{iR+B,iC+B}]; disp(['stitch: ',fName]);
            fData = load([directory_s,'\',fName,'.mat']);
            clear mask;
            if (iE~=iE_ref)
                try
                    mask = logical(fData.mask);
                    fData.sigma(~mask) = -1;
                    disp('have mask');
                catch
                end
            end
            
            switch VarLetter
                case 'xyuv'             % this works on the raw DIC data
                    u_conv = fData.u;
                    v_conv = fData.v;
                    exx_conv = fData.exx;
                    exy_conv = fData.exy;
                    eyy_conv = fData.eyy;
                case 'XYUV'             % this works best for the data after we performed distortion correction
                    u_conv = fData.U;
                    v_conv = fData.V;
                    exx_conv = fData.exx;
                    exy_conv = fData.exy;
                    eyy_conv = fData.eyy;
                    fData.x = fData.X;
                    fData.y = fData.Y;
                case 'Lagrange'         % this is the data after we perform strain calculation,
                    u_conv = fData.U;
                    v_conv = fData.V;
                    exx_conv = fData.exx_Lagrange;
                    exy_conv = fData.exy_Lagrange;
                    eyy_conv = fData.eyy_Lagrange;
                    fData.x = fData.X;
                    fData.y = fData.Y;
            end
            
            u_conv(fData.sigma==-1)=NaN;
            v_conv(fData.sigma==-1)=NaN;
            exx_conv(fData.sigma==-1)=NaN;
            exy_conv(fData.sigma==-1)=NaN;
            eyy_conv(fData.sigma==-1)=NaN;
            
            if USE_AVG_FILTER
                avgFilter = ones(Filter_Size)/Filter_Size/Filter_Size;
                u_conv = conv2(fData.u, avgFilter, 'same');
                v_conv = conv2(fData.v, avgFilter, 'same');
                exx_conv = conv2(fData.exx, avgFilter, 'same');
                exy_conv = conv2(fData.exy, avgFilter, 'same');
                eyy_conv = conv2(fData.eyy, avgFilter, 'same');
            end
            
            % here use row/col_begin, instead of RS(1) and CS(1), because its overall begin 
            if CUTEDGE
                if (iR>row_begin)&&(iC>col_begin)
                    fData.sigma(1:size(exx_conv,1)-round((transY(iR+B,iC+B)-transY(iR+B-1,iC+B)+Oly)/dic_step), :) = -1;
                    fData.sigma(:, 1:size(exx_conv,2)-round((transX(iR+B,iC+B)-transX(iR+B,iC+B-1)+Oly)/dic_step)) = -1;
                    u_conv(1:size(exx_conv,1)-round((transY(iR+B,iC+B)-transY(iR+B-1,iC+B)+Oly)/dic_step), :) = NaN;
                    u_conv(:, 1:size(exx_conv,2)-round((transX(iR+B,iC+B)-transX(iR+B,iC+B-1)+Oly)/dic_step)) = NaN;
                    v_conv(1:size(exx_conv,1)-round((transY(iR+B,iC+B)-transY(iR+B-1,iC+B)+Oly)/dic_step), :) = NaN;
                    v_conv(:, 1:size(exx_conv,2)-round((transX(iR+B,iC+B)-transX(iR+B,iC+B-1)+Oly)/dic_step)) = NaN;
                elseif (iC>col_begin)     % fovR = first row, cut some on the left
                    fData.sigma(:, 1:size(exx_conv,2)-round((transX(iR+B,iC+B)-transX(iR+B,iC+B-1)+Oly)/dic_step)) = -1;
                    u_conv(:, 1:size(exx_conv,2)-round((transX(iR+B,iC+B)-transX(iR+B,iC+B-1)+Oly)/dic_step)) = NaN;
                    v_conv(:, 1:size(exx_conv,2)-round((transX(iR+B,iC+B)-transX(iR+B,iC+B-1)+Oly)/dic_step)) = NaN;
                elseif (iR>row_begin)     % fovC = first col, cut some on the top
                    fData.sigma(1:size(exx_conv,1)-round((transY(iR+B,iC+B)-transY(iR+B-1,iC+B)+Oly)/dic_step), :) = -1;
                    u_conv(1:size(exx_conv,1)-round((transY(iR+B,iC+B)-transY(iR+B-1,iC+B)+Oly)/dic_step), :) = NaN;
                    v_conv(1:size(exx_conv,1)-round((transY(iR+B,iC+B)-transY(iR+B-1,iC+B)+Oly)/dic_step), :) = NaN;
                end
            end
            
            % find index of meaningful data points in global and local coord
            ind_patch_local = (fData.sigma~=-1);    % Locally meaningful data (simga~=-1)
            [nR,nC] = size(ind_patch_local);
            
            % YOU have one fov(r1c10) of different dimension! So, I need to round the data. Forget which specific test+FOV this is, but keep just in case.
            %             fData.x = dic_step * floor(fData.x/dic_step);
            %             fData.y = dic_step * floor(fData.y/dic_step);
            
            Ci = find(X(1,:) == fData.x(1)+transX(iR+B,iC+B));   % find starting point of patch in whole matrix.
            Ri = find(Y(:,1) == fData.y(1)+transY(iR+B,iC+B));
            Cf = Ci+nC-1;
            Rf = Ri+nR-1;
            ind_patch_global = zeros(size(X));
            ind_patch_global(Ri:Rf,Ci:Cf)=ind_patch_local;    % globally meaningful data (sigma ~=-1)
            
            ind_patch_global = logical(ind_patch_global);     % convert to logical
            ind_patch_local = logical(ind_patch_local);
            
            % sigma
            sigmaPatch = ones(size(X))*nan;     % initialize
            sigmaPatch(ind_patch_global) = fData.sigma(ind_patch_local);    % copy to global position
            if (iR==RS(1))&&(iC==CS(1))
                sigma = sigmaPatch;             % special case, 1st fov
            end
            
            ind_intercept_FOV = ~isnan(sigma) & ~isnan(sigmaPatch);     % ind_intercept labels the intercept region of the valid data.
            
            % for global method, we already determined how to adjust u/v. The reason to find intercept region is that, we don't want to cut too much to leave holes.
            % for local method, we need overlap regions with correlated u/v data, i.e., sigma~=-1, to determine the adjustment.  
            if strcmpi(stitch_method,'global')
                ind_intercept = ind_intercept_FOV;
            elseif strcmpi(stitch_method,'local')
                % This seems like a historic issue, but not solved yet. And maybe, actually, this doesn't need to be solved. 
                % Right now, this operation does not make a difference to ind_intercept.
                % Potentially, some more data can be labeled not to be used, so that ind_intercept labels the intercept region to be used to calculate uv adjustment.
                ind_intercept = ~isnan(sigma) & (sigmaPatch>-1);       
            end
        end
        
        ind_strain = ind_patch_global & (~ind_intercept_FOV);       % meaningful data, and not overlapped data. This is the region where new strain data to be pasted onto.
        if(iR==RS(1))&&(iC==CS(1))      
            ind_strain = ind_patch_global;      % special case, 1st image
        end
        ind_strain_local = ind_strain(Ri:Rf,Ci:Cf);
        
        sigma = nanmean(cat(3,sigma,sigmaPatch),3);
        clear sigmaPatch;
        
        % Considering the historic issue, a summary: 
        % ind_patch = ind_intercept_FOV + ind_strain.
        % ind_intercept_FOV = ind_intercept (stitched globally)
        % ind_intercept_FOV = ind_intercept + addtional regions to throw away (stitched locally) 
        
        % U
        uPatch = zeros(size(X))*nan;
        uPatch(ind_patch_global) = u_conv(ind_patch_local);
        if(iR==RS(1))&&(iC==CS(1))
            u = uPatch;
        end
        
        if strcmpi(stitch_method,'local')
            u1 = u(ind_intercept);
            u1 = nanmean(u1(:));
            u2 = uPatch(ind_intercept);
            u2 = nanmean(u2(:));
            uTrans(iR+B,iC+B) = u1 - u2;            
        end
        
        uPatch = uPatch + uTrans(iR+B,iC+B);
        u = nanmean(cat(3,u,uPatch),3);
        clear uPatch;
        
        % V
        vPatch = zeros(size(X))*nan;
        vPatch(ind_patch_global) = v_conv(ind_patch_local);
        if (iR==RS(1))&&(iC==CS(1))
            v = vPatch;
        end
        
        if strcmpi(stitch_method,'local')
            v1 = v(ind_intercept);
            v1 = nanmean(v1(:));
            v2 = vPatch(ind_intercept);
            v2 = nanmean(v2(:));
            vTrans(iR+B,iC+B) = v1 - v2;
        end
        
        vPatch = vPatch + vTrans(iR+B,iC+B);
        v = nanmean(cat(3,v,vPatch),3);
        clear vPatch;
        
        % exx, exy, eyy, simply fill with the new one
        exx(ind_strain) = exx_conv(ind_strain_local);
        exy(ind_strain) = exy_conv(ind_strain_local);
        eyy(ind_strain) = eyy_conv(ind_strain_local);
    end
    sigma(isnan(sigma)) = -1;
    
    if strcmpi(stitch_method,'local')
        save([directory_n,'\','methodLocal_uvTrans_e',STOP{iE+B}],'uTrans','vTrans');
    end
    
    if reduce_step ~= 1
        % data can be taken in reduced density
        x=X(1:reduce_step:end,1:reduce_step:end);
        y=Y(1:reduce_step:end,1:reduce_step:end);
        u=u(1:reduce_step:end,1:reduce_step:end);
        v=v(1:reduce_step:end,1:reduce_step:end);
        exx=exx(1:reduce_step:end,1:reduce_step:end);
        exy=exy(1:reduce_step:end,1:reduce_step:end);
        eyy=eyy(1:reduce_step:end,1:reduce_step:end);
        sigma=sigma(1:reduce_step:end,1:reduce_step:end);
    else
        x=X;
        y=Y;
    end
    
    if strcmpi(stitch_method,'global')
        save([directory_n,'\','global_method_stitched',f2,STOP{iE+B}],'x','y','u','v','exx','exy','eyy','sigma');
    elseif strcmpi(stitch_method,'local')
        save([directory_n,'\','local_method_stitched',f2,STOP{iE+B}],'x','y','u','v','exx','exy','eyy','sigma');
    end
    disp([directory_n,'\',f2,STOP{iE+B},' done']);
    
end

