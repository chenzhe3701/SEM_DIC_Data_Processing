% interpolate data using 'fit' or 'interp2' method
% Dt = interp_data(xf,yf,Df,xt,yt,tform_fwd,method, interpMethod)
%
% chenzhe, 2017-05-31
% chenzhe, 2018-05-25, update to solve possible problems with EBSD
% wrapping unwrapping. Use with caution.

function Dt = interp_data(xf,yf,Df,xt,yt,tform_fwd,method,interpMethod)
% if empty, no transform, just interp

% correct EBSD unwrapping related issue. Use with caution.
for iR = size(xf,1):-1:1
    if ~isequal(xf(iR,:),xf(1,:),'rows')
        xf(iR,:) = xf(1,:);
    end
end
yStep = yf(2)-yf(1);
for iR = 3:size(yf,1)
   if length(unique(yf(iR,:))>1)
      yf(iR,:) = yf(iR-1,:)+yStep; 
   end
end
    
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