% append deformed data by interpolation.
% Vic2D has a variable to select according to its documentation.
% Xp = X + U.  But this doesn't seem to work.
%
% Here, I define Xp as (X+U) interpolated to grid data.
%
% chenzhe, 2017-06-09, change the eij interp method from 'linear' to 'nearest'
%
% chenzhe, 2017-09-20, add note and edit
% (1) Not all data will have misalignment, so correct rotation-corrected
% data is not necessary
% (2) It's recommended to have a mask, otherwise all region including
% regions outside of the sample can be interpolated.
% To make a mask, consider other codes in this folder. 
% However, if there is not a mask available, this code just create one
% using sigma to let the code run.
% (3) By default, correct data with regular names, including x, y, u, v,
% exx, exy, eyy, because not all data should have rotation correction
% problem.
% (4) But if you have rotation corrected data, such as xR, vR, exyR, also
% correct it.
% (5) changed the variable names 'Xp'->'xp', 'Yp'->'yp'

[f,p] = uigetfile('D:\','Select .mat files to process','multiselect','on');
if ~iscell(f)
    f = cellstr(f);
end
[fm,pm] = uigetfile('D:\','Select the mask for non-rotation-corrected data. Recommend.  If do not have, just click cancel.');
[fmR,pmR] = uigetfile('D:\','Select the mask for rotation-corrected data. Recommend.  If do not have, just click cancel.');

%%
for ii = 1:length(f)
    % (1) by defaut, work on these variable
    data = load([p,f{ii}]);
    
    x = data.x;
    y = data.y;
    u = data.u;
    v = data.v;
    exx = data.exx;
    eyy = data.eyy;
    exy = data.exy;
    sigma = data.sigma;
    
    clear data;
    
    stepSize = x(1,2) - x(1,1);
    xd = x+u;
    yd = y+v;
    ad_x = mod(min(xd(:)) - min(x(:)), stepSize);   % some adjustment
    ad_y = mod(min(yd(:)) - min(y(:)), stepSize);
    a = floor(min(xd(:)) - ad_x);
    b = ceil(max(xd(:)) + ad_x);
    c = floor(min(yd(:)) - ad_y);
    d = ceil(max(yd(:)) + ad_y);
    [xp,yp] = meshgrid(a:stepSize:b, c:stepSize:d); % a meshgrid for deformed position
    [nR,nC] = size(xp);
    ind = (sigma(:)>-1);
    save([p,f{ii}],'xp','yp','-append');
    
    % This method is slow
    %     Ffit = fit([xp(ind),yp(ind)],exx(ind),'linearinterp');
    %     exxFit = Ffit(xp,yp);
    
    % make 'simgap' very similar to the original 'sigma'.  Then, based on
    % this sigmap, correct the interpolated data later, i.e., cut the
    % regions outside of the sample.
    % --> However, this turned out not good. sigmap can also go outside of
    % the sample region.
    %     F = scatteredInterpolant(xd(ind),yd(ind),sigma(ind),'nearest','none');
    %     sigmap = F(xp,yp);
    %     save([p,f{ii}],'sigmap','-append');
    %     clear sigmap;
    
    % First, interpolate u/v, to interpolate a sigma, that defines the
    % valid sample region. The purpose is just to get this sigma.
    % Then, use the valid u/v to interpolate.  It's also possible to
    % interpolate using the interpolated u/v (i.e., uu/vv).  To do this,
    % just make Xd = Xdd, Yd = Ydd.
    % However, possibly some u/v are interpolated bad, so the sigma still
    % looks not so good.
    %     F = scatteredInterpolant(x(ind),y(ind),u(ind),'linear');
    %     uu = F(x,y);
    %     F = scatteredInterpolant(x(ind),y(ind),v(ind),'linear');
    %     vv = F(x,y);
    %     Xdd = x+uu;
    %     Ydd = y+vv;
    %     F = scatteredInterpolant(Xdd(:),Ydd(:),sigma(:),'nearest','none');
    %     sigmap = F(xp,yp);
    %     UseInterpolatedUV = 0;
    %     if UseInterpolatedUV==1
    %         xd = Xdd;
    %         yd = Ydd;
    %     else
    %     end
    %
    %     clear uu vv Xdd Ydd;
    
    % So, the solution is to ALSO use the mask.  Interpolate the mask
    try
        mask = load([pm,fm]);
        mask = mask.mask;        %mask = mask.mask_R;
        indP = (mask>0);
        
        maskp = zeros(size(xp));
        cc = round((xd(indP)-xp(1))/stepSize);
        rr = round((yd(indP)-yp(1))/stepSize);
        for kk = 1:length(cc)
            try
                maskp(rr(kk),cc(kk))=1;
            catch
            end
        end
        
        maskp = grow_boundary(maskp);
        maskp = shrink_boundary(maskp);
        
%         % this is the old method, but I found it sometimes doesn't  work well on 2017-09-20, so I replaced it with a new/maybe also faster method. 
%         F = scatteredInterpolant(xd(ind),yd(ind),mask(ind),'linear','none');
%         maskp = F(xp,yp);
%         maskp(maskp>0.999999)=1;
%         maskp(maskp~=1)=0;
    catch
        disp('no mask selected, creating one using sigma value');
        % create it using xd,yd, sigma, xp, yp
        mask = sigma;
        
        maskp = zeros(size(xp));
        cc = round((xd(ind)-xp(1))/stepSize);
        rr = round((yd(ind)-yp(1))/stepSize);
        for kk = 1:length(cc)
            try
                maskp(rr(kk),cc(kk))=1;
            catch
            end
        end
    
        maskp = grow_boundary(maskp);
        maskp = shrink_boundary(maskp);
    
    end
    
    try
        save([pm,fm],'maskp','-append');    % append the deformed mask
    catch
    end
    disp('interpolated mask');
    
    F = scatteredInterpolant(xd(ind),yd(ind),exx(ind),'nearest','none');
    exxp = F(xp,yp);
    exxp(maskp~=1) = nan;
    save([p,f{ii}],'exxp','-append');    % Note: F is too big to save.  Not worth to save.
    clear exxp;
    disp('interpolated exx');
    
    F = scatteredInterpolant(xd(ind),yd(ind),exy(ind),'nearest','none');
    exyp = F(xp,yp);
    exyp(maskp~=1) = nan;
    save([p,f{ii}],'exyp','-append');
    clear exyp;
    disp('interpolated exy');
    
    F = scatteredInterpolant(xd(ind),yd(ind),eyy(ind),'nearest','none');
    eyyp = F(xp,yp);
    eyyp(maskp~=1) = nan;
    save([p,f{ii}],'eyyp','-append');
    clear eyyp;
    disp('interpolated eyy');
    
    F = scatteredInterpolant(xd(ind),yd(ind),u(ind),'linear','none');
    up = F(xp,yp);
    up(maskp~=1) = nan;
    save([p,f{ii}],'up','-append');
    clear up;
    disp('interpolated u');
    
    F = scatteredInterpolant(xd(ind),yd(ind),v(ind),'linear','none');
    vp = F(xp,yp);
    vp(maskp~=1) = nan;
    save([p,f{ii}],'vp','-append');
    clear vp;
    disp('interpolated v');
    
    F = scatteredInterpolant(xd(ind),yd(ind),sigma(ind),'nearest','none');
    sigmap = F(xp,yp);
    sigmap(maskp~=1)=-1;
    save([p,f{ii}],'sigmap','-append');
    clear sigmap;
    disp('interpolated sigma');
    
    % (2) see if there is rotation corrected data to work on
    try
        data = load([p,f{ii}]);
        
        xR = data.xR;
        yR = data.yR;
        uR = data.uR;
        vR = data.vR;
        exxR = data.exxR;
        eyyR = data.eyyR;
        exyR = data.exyR;
        sigmaR = data.sigmaR;
        
        clear data;
        
        stepSize = xR(1,2) - xR(1,1);
        xRd = xR+uR;
        yRd = yR+vR;
        ad_x = mod(min(xRd(:)) - min(xR(:)), stepSize);   % some adjustment
        ad_y = mod(min(yRd(:)) - min(yR(:)), stepSize);
        a = floor(min(xRd(:)) - ad_x);
        b = ceil(max(xRd(:)) + ad_x);
        c = floor(min(yRd(:)) - ad_y);
        d = ceil(max(yRd(:)) + ad_y);
        [xRp,yRp] = meshgrid(a:stepSize:b, c:stepSize:d); % a meshgrid for deformed position
        [nR,nC] = size(xRp);
        ind = (sigmaR(:)>-1);
        save([p,f{ii}],'xRp','yRp','-append');
        
               
        try
            maskR = load([pmR,fmR]);
            maskR = maskR.maskR;
            indP = (maskR>0);
            
            maskRp = zeros(size(xRp));
            cc = round((xRd(indP)-xRp(1))/stepSize);
            rr = round((yRd(indP)-yRp(1))/stepSize);
            for kk = 1:length(cc)
                try
                    maskRp(rr(kk),cc(kk))=1;
                catch
                end
            end
            
            maskRp = grow_boundary(maskRp);
            maskRp = shrink_boundary(maskRp);
        
%             F = scatteredInterpolant(xRd(ind),yRd(ind),maskR(ind),'linear','none');
%             maskRp = F(xRp,yRp);
%             maskRp(maskRp>0.999999)=1;
%             maskRp(maskRp~=1)=0;
        catch
            disp('no mask selected, creating one using sigmaR value');
            % create it using xd,yd, sigma, xp, yp
            maskR = sigmaR;
            
            maskRp = zeros(size(xRp));
            cc = round((xRd(ind)-xRp(1))/stepSize);
            rr = round((yRd(ind)-yRp(1))/stepSize);
            for kk = 1:length(cc)
                try
                    maskRp(rr(kk),cc(kk))=1;
                catch
                end
            end
            
            maskRp = grow_boundary(maskRp);
            maskRp = shrink_boundary(maskRp);
            
        end

        try
            save([pmR,fmR],'maskRp','-append');    % append the deformed mask
        catch
        end
        disp('interpolated maskR');
        
        F = scatteredInterpolant(xRd(ind),yRd(ind),exxR(ind),'nearest','none');
        exxRp = F(xRp,yRp);
        exxRp(maskRp~=1) = nan;
        save([p,f{ii}],'exxRp','-append');    % Note: F is too big to save.  Not worth to save.
        clear exxRp;
        disp('interpolated exxR');
        
        F = scatteredInterpolant(xRd(ind),yRd(ind),exyR(ind),'nearest','none');
        exyRp = F(xRp,yRp);
        exyRp(maskRp~=1) = nan;
        save([p,f{ii}],'exyRp','-append');
        clear exyRp;
        disp('interpolated exyR');
        
        F = scatteredInterpolant(xRd(ind),yRd(ind),eyyR(ind),'nearest','none');
        eyyRp = F(xRp,yRp);
        eyyRp(maskRp~=1) = nan;
        save([p,f{ii}],'eyyRp','-append');
        clear eyyRp;
        disp('interpolated eyyR');
        
        F = scatteredInterpolant(xRd(ind),yRd(ind),uR(ind),'linear','none');
        uRp = F(xRp,yRp);
        uRp(maskRp~=1) = nan;
        save([p,f{ii}],'uRp','-append');
        clear uRp;
        disp('interpolated uR');
        
        F = scatteredInterpolant(xRd(ind),yRd(ind),vR(ind),'linear','none');
        vRp = F(xRp,yRp);
        vRp(maskRp~=1) = nan;
        save([p,f{ii}],'vRp','-append');
        clear vRp;
        disp('interpolated vR');
        
        F = scatteredInterpolant(xRd(ind),yRd(ind),sigmaR(ind),'nearest','none');
        sigmaRp = F(xRp,yRp);
        sigmaRp(maskRp~=1)=-1;
        save([p,f{ii}],'sigmaRp','-append');
        clear sigmaRp;
        disp('interpolated sigmaR');
        
    catch
        clear data;
    end
    
    disp(['finished file',f{ii}]);
end

