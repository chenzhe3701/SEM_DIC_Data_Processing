% append deformed data by interpolation.
% Vic2D has a variable to select according to its documentation.
% Xp = X + U.  But this doesn't seem to work.
%
% Here, I define Xp as (X+U) interpolated to grid data.
%
% chenzhe, 2017-06-09, change the eij interp method from 'linear' to 'nearest' 

[f,p] = uigetfile('D:\','Select .mat files','multiselect','on');
if ~iscell(f)
    f = cellstr(f);
end
[fm,pm] = uigetfile('D:\','Select the mask');

for ii = 1:length(f)
    data = load([p,f{ii}]);
    try
        x = data.xR;
        y = data.yR;
        u = data.uR;
        v = data.vR;
        exx = data.exxR;
        eyy = data.eyyR;
        exy = data.exyR;
        sigma = data.sigmaR;
    catch
        x = data.x;
        y = data.y;
        u = data.u;
        v = data.v;
        exx = data.exx;
        eyy = data.eyy;
        exy = data.exy;
        sigma = data.sigma;
    end
    clear data;
    
    stepSize = x(1,2) - x(1,1);
    Xd = x+u;
    Yd = y+v;
    ad_x = mod(min(Xd(:)) - min(x(:)), stepSize);
    ad_y = mod(min(Yd(:)) - min(y(:)), stepSize);
    a = floor(min(Xd(:)) - ad_x);
    b = ceil(max(Xd(:)) + ad_x);
    c = floor(min(Yd(:)) - ad_y);
    d = ceil(max(Yd(:)) + ad_y);
    [Xp,Yp] = meshgrid(a:stepSize:b, c:stepSize:d);
    [nR,nC] = size(Xp);
    ind = (sigma(:)>-1);
    save([p,f{ii}],'Xp','Yp','-append');
    
    % This method is slow
    %     Ffit = fit([Xp(ind),Yp(ind)],exx(ind),'linearinterp');
    %     exxFit = Ffit(Xp,Yp);
  
    % make 'simgap' very similar to the original 'sigma'.  Then, based on
    % this sigmap, correct the interpolated data later, i.e., cut the
    % regions outside of the sample.
    % --> However, this turned out not good. sigmap can also go outside of
    % the sample region.    
    %    F = scatteredInterpolant(Xd(ind),Yd(ind),sigma(ind),'nearest','none');
    %    sigmap = F(Xp,Yp);
    %    save([p,f{ii}],'sigmap','-append');
    %    clear sigmap;
    
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
    %     sigmap = F(Xp,Yp);
    %     UseInterpolatedUV = 0;
    %     if UseInterpolatedUV==1
    %         Xd = Xdd;
    %         Yd = Ydd;
    %     else
    %     end
    %
    % clear uu vv Xdd Ydd;
    
    % So, the solution is to ALSO use the mask.  Interpolate the mask
    mask = load([pm,fm]);
    try 
        mask = mask.mask;
    catch
        mask = mask.mask_R;
    end
    F = scatteredInterpolant(Xd(ind),Yd(ind),mask(ind),'linear','none');
    maskp = F(Xp,Yp);
    maskp(maskp>0.999999)=1;
    maskp(maskp~=1)=0;
    save([pm,fm],'maskp','-append');    % append the deformed mask
    disp('interpolated mask');
    
    F = scatteredInterpolant(Xd(ind),Yd(ind),exx(ind),'nearest','none');
    exxp = F(Xp,Yp);
    exxp(maskp~=1) = nan;
    save([p,f{ii}],'exxp','-append');    % Note: F is too big to save.  Not worth to save.
    clear exxp;
    disp('interpolated exx');
    
    F = scatteredInterpolant(Xd(ind),Yd(ind),exy(ind),'nearest','none');
    exyp = F(Xp,Yp);
    exyp(maskp~=1) = nan;
    save([p,f{ii}],'exyp','-append');
    clear exyp;
    disp('interpolated exy');
    
    F = scatteredInterpolant(Xd(ind),Yd(ind),eyy(ind),'nearest','none');
    eyyp = F(Xp,Yp);
    eyyp(maskp~=1) = nan;
    save([p,f{ii}],'eyyp','-append');
    clear eyyp;
    disp('interpolated eyy');
    
    F = scatteredInterpolant(Xd(ind),Yd(ind),u(ind),'linear','none');
    up = F(Xp,Yp);
    up(maskp~=1) = nan;
    save([p,f{ii}],'up','-append');
    clear up;
    disp('interpolated u');
    
    F = scatteredInterpolant(Xd(ind),Yd(ind),v(ind),'linear','none');
    vp = F(Xp,Yp);
    vp(maskp~=1) = nan;
    save([p,f{ii}],'vp','-append');
    clear vp;
    disp('interpolated v');
    
    F = scatteredInterpolant(Xd(ind),Yd(ind),sigma(ind),'nearest','none');
    sigmap = F(Xp,Yp);
    sigmap(maskp~=1)=-1;
    save([p,f{ii}],'sigmap','-append');
    clear sigmap;
    disp('interpolated sigma');
    
    disp(['finished file',f{ii}]);
end

