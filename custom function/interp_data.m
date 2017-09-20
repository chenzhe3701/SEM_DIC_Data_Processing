% interpolate data using 'fit' or 'interp2' method
% Dt = interp_data(xf,yf,Df,xt,yt,tform_fwd,method, interpMethod)
%
% chenzhe, 2017-05-31

function Dt = interp_data(xf,yf,Df,xt,yt,tform_fwd,method,interpMethod)
% if empty, no transform, just interp
if isempty(tform_fwd)
    tform_fwd = maketform('projective',diag([1 1 1]));
end
switch method
    case {'interp','Interp'}
        [xt_bwd, yt_bwd] = tforminv(tform_fwd,xt,yt);
        Dt = interp2(xf,yf,Df,xt_bwd,yt_bwd,interpMethod);
    case {'fit','Fit'}
        [xf_fwd, yf_fwd] = tformfwd(tform_fwd,[xf(:),yf(:)]);
        F = scatteredInterpolant(xf_fwd, yf_fwd, Df(:),interpMethod);
        Dt = F(xt,yt);
end
end