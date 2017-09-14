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

[f,p] = uigetfile('D:\Marissa_test_20170430_renamed_cropped_sequential_rotated\translations_searched_vertical_stop_0.mat','select translation');
load([p,f]);    % load translation data
REFSTOP = 1;    % stop number for reference images, default=1

dic_step = 5;       % step size in DIC analysis
try
    % x0 = xTransFOV
    transX = transX
catch
    transX = [];    % can manual input
end
try
    % y0 = yTransFOV
    transY = transY
catch
    transY = [];
end
transX = dic_step * round(transX/dic_step);
transY = dic_step * round(transY/dic_step);
transX_incremental = dic_step * round(transX_incremental/dic_step);
transY_incremental = dic_step * round(transY_incremental/dic_step);


% FOV = {
%     '001'   '002'   '003'   '004'   '005'   '006'
%     '012'   '011'   '010'   '009'   '008'   '007'
%     '013'	'014'	'015'	'016'	'017'	'018'
%     '024'	'023'	'022'	'021'	'020'	'019'
%     '025'	'026'	'027'	'028'	'029'	'030'
%     '036'	'035'	'034'	'033'	'032'	'031'};
% Marrisa you may use this?
% FOV = {
% 'r0c0', 'r0c1', 'r0c2', 'r0c3', 'r0c4', 'r0c5', 'r0c6', 'r0c7', 'r0c8', 'r0c9', 'r0c10', 'r0c11','r0c12', 'r0c13';
% 'r1c0', 'r1c1', 'r1c2', 'r1c3', 'r1c4', 'r1c5', 'r1c6', 'r1c7', 'r1c8', 'r1c9', 'r1c10', 'r1c11','r1c12', 'r1c13';
% 'r2c0', 'r2c1', 'r2c2', 'r2c3', 'r2c4', 'r2c5', 'r2c6', 'r2c7', 'r2c8', 'r2c9', 'r2c10', 'r2c11','r2c12', 'r2c13';
% 'r3c0', 'r3c1', 'r3c2', 'r3c3', 'r3c4', 'r3c5', 'r3c6', 'r3c7', 'r3c8', 'r3c9', 'r3c10', 'r3c11','r3c12', 'r3c13';
% };
minRNum = 0; % starting # of FOV rows
maxRNum = 3;
minCNum = 0;
maxCNum = 13;  % ending # of FOV cols
for iR = minRNum:maxRNum
    for iC = minCNum:maxCNum
        FOV{iR+minRNum+1,iC+minCNum+1} = ['r',num2str(iR),'c',num2str(iC)];
    end
end
nR_FOV = maxRNum-minRNum+1;
nC_FOV = maxCNum-minCNum+1;

% f1 = 'T5_#7_fov_'; f2='_stop_'; f3='.mat';
% Marissa you may use this?
f1 = '20170430_ts5Al_02_e'; f2='_';

% STOP = {'001','002','003','004','005','006','007'};
% STOP = {'001','002','003','004','005','006','007','008','009','010','011','012','013','014','015'};
% Marissa you may try+change this
STOP = {'0','1','2','3','4','5','6'};

directory_s = uigetdir('D:\Marissa_test_20170430\4]_20170430_ts5Al_02 tensile test_non_incremental_DIC_rotated','choose a single directory that directly contains all the DIC mat files');
directory_n = uigetdir('D:\','choose the new/destination directory');

% dic_size = 21800;   % total width of all images.
% [X,Y] = meshgrid(1:dic_step:dic_size);    % you can design the size of the data matrix

% Marissa you may want to try the following instead ------------------------- NOTE NOTE NOTE
ii = 0;     % basically, this is x(1) in each FOV.
resX = 6144; 
resY = 4096;
dic_size_x = 70000;
dic_size_y = 14500;
dataStartX = ii + min(transX(:));   dic_size_x = max(dic_size_x, dataStartX + max(transX(:))+resX)
dataStartY = ii + min(transY(:));   dic_size_y = max(dic_size_y, dataStartY + max(transY(:))+resY)
[X,Y] = meshgrid(dataStartX:dic_step:dic_size_x, dataStartY:dic_step:dic_size_y);

fovR_start = 1;     % range in the 'FOV' that you want to merge.
fovR_stop = 4;
fovC_start = 1;
fovC_stop = 14;

Filter_Size = 1;    % odd number, smooth with avg filter, then reduce matrix size.  This can smooth DIC data.
reduce_step = 1;    % reduce data size, data point step size.
%%%%%%

% use avg_filter, in order to reduce size later
USE_AVG_FILTER = 0;
% 'xyuv' for raw DIC data, 'XYUV' for distortion corrected, 'Lagrange' for
% distortion corrected and strain calculated.
VarLetter = 'xyuv';
CUTEDGE = 1;        % cut edge or not?
OVERLAY = 200;      % how many overlay to give
nSamples = 10;      % how many points in each row/column to determine u/v shift
%%
for iStop = 1:length(STOP)
    exx = zeros(size(X))*nan;
    exy = zeros(size(X))*nan;
    eyy = zeros(size(X))*nan;
    sigma = ones(size(X))*nan;
    u = zeros(size(X))*nan;
    v = zeros(size(X))*nan;
    
    N = nR_FOV*nC_FOV;  % total number of elements/unknown translations/FOVs
    % Y = (K.*A)*X
    A = zeros(N^2,N);
    Yu = zeros(N^2,1);
    Yv = zeros(N^2,1);
    K = zeros(N^2,1);
    
    B1 = zeros(N,N);
    B2 = zeros(1,N);
    B3 = 0;
    D1 = zeros(N,N);
    D2 = zeros(1,N);
    D3 = 0;
    
    for ii = 1:N                % ii = 1st fov #
        for jj = ii+1:N         % jj = 2nd fov #
            rr = (ii-1)*N + jj;
            A(rr,ii) = 1;
            A(rr,jj) = -1;
        end
    end
    
    for fovR = fovR_start:fovR_stop
        for fovC = fovC_start:fovC_stop
            % (0) current fov data
            fName = [f1,STOP{iStop},f2,FOV{fovR,fovC}]; disp(fName);
            fData = load([directory_s,'\',fName,'.mat']);
            clear mask;
            if (iStop~=REFSTOP)
                try
                    mask = logical(fData.mask);
                    fData.sigma(~mask) = -1;
                    disp('have mask');
                catch
                end
            end
            flag = double((fData.sigma~=-1)&(~isnan(fData.sigma)));    % use this to label overlap region by looking at sigma values
            
            [nR,nC] = size(fData.sigma);
            [~,Ci] = min(abs(fData.x(1)+transX(fovR,fovC) - X(1,:)));   % find starting point of patch in whole matrix.
            [~,Ri] = min(abs(fData.y(1)+transY(fovR,fovC) - Y(:,1)));
            Cf = Ci + nC - 1;
            Rf = Ri + nR - 1;
            
            % (1) search for right neighbor
            if fovC<fovC_stop
                fName_1 = [f1,STOP{iStop},f2,FOV{fovR,fovC+1}]; disp(fName_1);
                fData_1 = load([directory_s,'\',fName_1,'.mat']);
                clear mask;
                if (iStop~=REFSTOP)
                    try
                        mask = logical(fData_1.mask);
                        fData_1.sigma(~mask) = -1;
                        disp('have mask');
                    catch
                    end
                end
                flag_1 = double((fData_1.sigma~=-1)&(~isnan(fData_1.sigma)));    % use this to label overlap region by looking at sigma values
                
                [nR_1,nC_1] = size(fData_1.sigma);
                [~,Ci_1] = min(abs(fData_1.x(1)+transX(fovR,fovC+1) - X(1,:)));   % find starting point of patch in whole matrix.
                [~,Ri_1] = min(abs(fData_1.y(1)+transY(fovR,fovC+1) - Y(:,1)));
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
                
                flagT = flag;
                flagT(ri:rf,ci:cf) = flagT(ri:rf,ci:cf) + flag_1(ri_1:rf_1,ci_1:cf_1);
                flagT = double(flagT==2);   % This is the overlay region
                iii = find(sum(flagT,2),11);    % do not consider the 1st or last 10 rows/columns
                iii = iii(end);
                iij = find(sum(flagT,2),11,'last');
                iij = iij(1);
                for iRR = iii:iij
                    indC = find(flagT(iRR,:),nSamples,'last');
                    if ~isempty(indC)
                        indC = indC(1:min(nSamples,length(indC)));
                        flagT(iRR, indC) = 2;
                    end
                end
                flagT = (flagT==2); % logical
                flag_1 = false(size(flag_1));
                flag_1(ri_1:rf_1,ci_1:cf_1) = flagT(ri:rf,ci:cf);
                uShift = fData.u(flagT);
                uShift_1 = fData_1.u(flag_1);
                vShift = fData.v(flagT);
                vShift_1 = fData_1.v(flag_1);
                
                sub_1 = (fovR-1)*nC_FOV + fovC;
                sub_2 = (fovR-1)*nC_FOV + fovC+1;
                rr = (sub_1-1)*N + sub_2;
                
                K(rr) = K(rr)+length(uShift);
                Yu(rr) = Yu(rr) + sum(uShift_1) - sum(uShift);
                Yv(rr) = Yv(rr) + sum(vShift_1) - sum(vShift);
                
                B1(sub_1,sub_1) = B1(sub_1,sub_1) + length(uShift);
                B1(sub_2,sub_2) = B1(sub_2,sub_2) + length(uShift);
                B1(sub_1,sub_2) = B1(sub_1,sub_2) - 2*length(uShift);
                B2(sub_1) = B2(sub_1) + 2 * sum(uShift-uShift_1);
                B2(sub_2) = B2(sub_2) - 2 * sum(uShift-uShift_1);
                B3 = B3 + (uShift - uShift_1)'*(uShift - uShift_1);
                
                D1(sub_1,sub_1) = D1(sub_1,sub_1) + length(vShift);
                D1(sub_2,sub_2) = D1(sub_2,sub_2) + length(vShift);
                D1(sub_1,sub_2) = D1(sub_1,sub_2) - 2*length(vShift);
                D2(sub_1) = D2(sub_1) + 2 * sum(vShift-vShift_1);
                D2(sub_2) = D2(sub_2) - 2 * sum(vShift-vShift_1);
                D3 = D3 + (vShift - vShift_1)'*(vShift - vShift_1);
            end
            
            % (2) search for the bottom neighbor
            if fovR < fovR_stop
                fName_2 = [f1,STOP{iStop},f2,FOV{fovR+1,fovC}]; disp(fName_2);
                fData_2 = load([directory_s,'\',fName_2,'.mat']);
                clear mask;
                if (iStop~=REFSTOP)
                    try
                        mask = logical(fData_2.mask);
                        fData_2.sigma(~mask) = -1;
                        disp('have mask');
                    catch
                    end
                end
                flag_2 = double((fData_2.sigma~=-1)&(~isnan(fData_2.sigma)));    % use this to label overlap region by looking at sigma values
                
                [nR_2,nC_2] = size(fData_2.sigma);
                [~,Ci_2] = min(abs(fData_2.x(1)+transX(fovR+1,fovC) - X(1,:)));   % find starting point of patch in whole matrix.
                [~,Ri_2] = min(abs(fData_2.y(1)+transY(fovR+1,fovC) - Y(:,1)));
                Cf_2 = Ci_2 + nC_2 - 1;
                Rf_2 = Ri_2 + nR_2 - 1;
                
                ci = max(1, 1+Ci_2-Ci);
                cf = min(nC, nC + Cf_2-Cf);
                ri = max(1, 1+Ri_2-Ri);
                rf = min(nR, nR + Rf_2-Rf);
                
                ci_2 = ci + Ci-Ci_2;
                cf_2 = cf + Ci-Ci_2;
                ri_2 = ri + Ri-Ri_2;
                rf_2 = rf + Ri-Ri_2;
                
                flagT = flag;
                flagT(ri:rf,ci:cf) = flagT(ri:rf,ci:cf) + flag_2(ri_2:rf_2,ci_2:cf_2);
                flagT = double(flagT==2);   % This is the overlay region
                iii = find(sum(flagT,1),301);   % throw away left data
                iii = iii(end);
                iij = find(sum(flagT,1),11,'last');
                iij = iij(1);
                for iCC = iii:iij
                    indR = find(flagT(:,iCC),nSamples,'last');
                    if ~isempty(indR)
                        indR = indR(1:min(nSamples,length(indR)));
                        flagT(indR,iCC) = 2;
                    end
                end
                flagT = (flagT==2); % logical
                flag_2 = false(size(flag_2));
                flag_2(ri_2:rf_2,ci_2:cf_2) = flagT(ri:rf,ci:cf);
                uShift = fData.u(flagT);
                uShift_2 = fData_2.u(flag_2);
                vShift = fData.v(flagT);
                vShift_2 = fData_2.v(flag_2);
                
                sub_1 = (fovR-1)*nC_FOV + fovC;
                sub_2 = (fovR-1+1)*nC_FOV + fovC;
                rr = (sub_1-1)*N + sub_2;
                
                K(rr) = K(rr)+length(uShift);
                Yu(rr) = Yu(rr) + sum(uShift_2) - sum(uShift);
                Yv(rr) = Yv(rr) + sum(vShift_2) - sum(vShift);
                
                B1(sub_1,sub_1) = B1(sub_1,sub_1) + length(uShift);
                B1(sub_2,sub_2) = B1(sub_2,sub_2) + length(uShift);
                B1(sub_1,sub_2) = B1(sub_1,sub_2) - 2*length(uShift);
                B2(sub_1) = B2(sub_1) + 2 * sum(uShift-uShift_2);
                B2(sub_2) = B2(sub_2) - 2 * sum(uShift-uShift_2);
                B3 = B3 + (uShift - uShift_2)'*(uShift - uShift_2);
                
                D1(sub_1,sub_1) = D1(sub_1,sub_1) + length(vShift);
                D1(sub_2,sub_2) = D1(sub_2,sub_2) + length(vShift);
                D1(sub_1,sub_2) = D1(sub_1,sub_2) - 2*length(vShift);
                D2(sub_1) = D2(sub_1) + 2 * sum(vShift-vShift_2);
                D2(sub_2) = D2(sub_2) - 2 * sum(vShift-vShift_2);
                D3 = D3 + (vShift - vShift_2)'*(vShift - vShift_2);                
            end
            disp(['RC: ',num2str(fovR),' , ',num2str(fovC)]);
        end
    end
    
    ind = (K~=0);
    K = K(ind);
    A = A(ind,:);
    Yu = Yu(ind);
    Yv = Yv(ind);
    KA = K.*A;
    
    uTrans = KA\Yu;
    uTrans = transpose(reshape(uTrans-uTrans(1),nC_FOV,nR_FOV));
    vTrans = KA\Yv;
    vTrans = transpose(reshape(vTrans-vTrans(1),nC_FOV,nR_FOV));
    
    uTrans1 = fminunc(@(x) sum(sum(B1.*(x*x')))+B2*x+B3, reshape(uTrans',[],1));
    uTrans1 = transpose(reshape(uTrans1-uTrans1(1),nC_FOV,nR_FOV));
    vTrans1 = fminunc(@(x) sum(sum(D1.*(x*x')))+D2*x+D3, reshape(uTrans',[],1));
    vTrans1 = transpose(reshape(vTrans1-vTrans1(1),nC_FOV,nR_FOV));
    
    % This doesn't work good
    uTrans2 = inv(KA'*KA)*KA'*Yu;
    uTrans2 = transpose(reshape(uTrans2-uTrans2(1),nC_FOV,nR_FOV));
    vTrans2 = inv(KA'*KA)*KA'*Yv;
    vTrans2 = transpose(reshape(vTrans2-vTrans2(1),nC_FOV,nR_FOV));
    save([directory_n,'\','uvTrans_e',STOP{iStop}],'uTrans','vTrans','uTrans1','vTrans1','uTrans2','vTrans2');
    
    for fovR = fovR_start:fovR_stop
        for fovC = fovC_start:fovC_stop
            ind_intercept = 0;
            Oly = OVERLAY-150 + 300 * iStop/length(STOP);
            while(sum(ind_intercept(:))==0)
                Oly = Oly + 150;
                % fName = [f1,FOV{fovR,fovC},f2,STOP{iStop},f3];
                % Marrisa try this?
                fName = [f1,STOP{iStop},f2,FOV{fovR,fovC}]; disp(['stitch: ',fName]);
                fData = load([directory_s,'\',fName,'.mat']);
                
                clear mask;
                if (iStop~=REFSTOP)
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
                
                if CUTEDGE
                    if (fovR>fovR_start)&&(fovC>fovC_start)
                        fData.sigma(1:size(exx_conv,1)-round((transY(fovR,fovC)-transY(fovR-1,fovC)+Oly)/dic_step), :) = -1;
                        fData.sigma(:, 1:size(exx_conv,2)-round((transX(fovR,fovC)-transX(fovR,fovC-1)+Oly)/dic_step)) = -1;
                        u_conv(1:size(exx_conv,1)-round((transY(fovR,fovC)-transY(fovR-1,fovC)+Oly)/dic_step), :) = NaN;
                        u_conv(:, 1:size(exx_conv,2)-round((transX(fovR,fovC)-transX(fovR,fovC-1)+Oly)/dic_step)) = NaN;
                        v_conv(1:size(exx_conv,1)-round((transY(fovR,fovC)-transY(fovR-1,fovC)+Oly)/dic_step), :) = NaN;
                        v_conv(:, 1:size(exx_conv,2)-round((transX(fovR,fovC)-transX(fovR,fovC-1)+Oly)/dic_step)) = NaN;
                    elseif (fovC>fovC_start)     % fovR = 1, cut some on the left
                        fData.sigma(:, 1:size(exx_conv,2)-round((transX(fovR,fovC)-transX(fovR,fovC-1)+Oly)/dic_step)) = -1;
                        u_conv(:, 1:size(exx_conv,2)-round((transX(fovR,fovC)-transX(fovR,fovC-1)+Oly)/dic_step)) = NaN;
                        v_conv(:, 1:size(exx_conv,2)-round((transX(fovR,fovC)-transX(fovR,fovC-1)+Oly)/dic_step)) = NaN;
                    elseif (fovR>fovR_start)     % fovC = 1
                        fData.sigma(1:size(exx_conv,1)-round((transY(fovR,fovC)-transY(fovR-1,fovC)+Oly)/dic_step), :) = -1;
                        u_conv(1:size(exx_conv,1)-round((transY(fovR,fovC)-transY(fovR-1,fovC)+Oly)/dic_step), :) = NaN;
                        v_conv(1:size(exx_conv,1)-round((transY(fovR,fovC)-transY(fovR-1,fovC)+Oly)/dic_step), :) = NaN;
                    end
                end
                
                % find index of meaningful data points in global and local coord
                ind_patch_local = (fData.sigma~=-1);    % Locally meaningful data (simga~=-1)
                [nR,nC] = size(ind_patch_local);
                
                % YOU have one fov(r1c10) of different dimension! So, I need to round the data
                %             fData.x = dic_step * floor(fData.x/dic_step);
                %             fData.y = dic_step * floor(fData.y/dic_step);
                
                C = find(X(1,:) == fData.x(1)+transX(fovR,fovC));   % find starting point of patch in whole matrix.  Note: this code allows some flexibility of the position (+-1 points, currently).
                R = find(Y(:,1) == fData.y(1)+transY(fovR,fovC));
                
                ind_patch_global = zeros(size(X));
                ind_patch_global(R:R+nR-1,C:C+nC-1)=ind_patch_local;    % globally meaningful data (sigma ~=-1)
                
                ind_patch_global = logical(ind_patch_global);     % convert to logical
                ind_patch_local = logical(ind_patch_local);
                
                % sigma
                sigmaPatch = ones(size(X))*nan;     % initialize
                sigmaPatch(ind_patch_global) = fData.sigma(ind_patch_local);    % copy to global position
                if (fovR==fovR_start)&&(fovC==fovC_start)
                    sigma = sigmaPatch;
                end
                
                ind_intercept = ~isnan(sigma) & ~isnan(sigmaPatch);     % ind_intercept labels the Intercept region of the valid data.
            end
            
            ind_strain = ind_patch_global & (~ind_intercept);       % meaningful data, and not overlapped data.
            if(fovR==fovR_start)&&(fovC==fovC_start)      % special case, 1st image
                ind_strain = ind_patch_global;
            end
            ind_strain_local = ind_strain(R:R+nR-1,C:C+nC-1);
            
            sigma = nanmean(cat(3,sigma,sigmaPatch),3);
            clear sigmaPatch;
            
            % we need to create a special intercept for u/v stitching, making it 'intercept' only on 1 edge, instead of 2 edges.
            % since we are not using [fData.sigma],[sigmaPatch] again, we just modify it, to reduce the overlap area
            
            
            % modify ind_intercept, to calculate U/V shift. Copy the 'local'  value, modify, then put it back.  Shrink the intercept region.
            % ind_intercept_local = ind_intercept(R:R+nR-1,C:C+nC-1);
            % ind_intercept(R:R+nR-1,C:C+nC-1) = ind_intercept_local;
            
            
            % U
            uPatch = zeros(size(X))*nan;
            uPatch(ind_patch_global) = u_conv(ind_patch_local);
            if(fovR==fovR_start)&&(fovC==fovC_start)
                u = uPatch;
            end
            
            uPatch = uPatch + uTrans(fovR,fovC);
            u = nanmean(cat(3,u,uPatch),3);
            clear uPatch;
            
            % V
            vPatch = zeros(size(X))*nan;
            vPatch(ind_patch_global) = v_conv(ind_patch_local);
            if (fovR==fovR_start)&&(fovC==fovC_start)
                v = vPatch;
            end
            vPatch = vPatch + vTrans(fovR,fovC);
            v = nanmean(cat(3,v,vPatch),3);
            clear vPatch;
            
            % exx, exy, eyy, simply fill with the new one
            exx(ind_strain) = exx_conv(ind_strain_local);
            exy(ind_strain) = exy_conv(ind_strain_local);
            eyy(ind_strain) = eyy_conv(ind_strain_local);
            
        end
    end
    sigma(isnan(sigma)) = -1;
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
    save([directory_n,'\',f2,STOP{iStop},'_oly_',num2str(OVERLAY)],'x','y','u','v','exx','exy','eyy','sigma');
    disp([directory_n,'\',f2,STOP{iStop},' done']);
    % Marissa you may try this
    % save([directory_n,'\',f1,STOP{iStop}],'x','y','u','v','exx','exy','eyy','sigma');
    
end

