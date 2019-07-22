% chenzhe, 2017-06-08.
% correct the exy field of DIC data, i.e., flip the sign.
% If corrected, append a varialbe exy_corrected = 1.
%
% chenzhe, 2017-09-19 note
% currently only corrects 'exy' variable.  Can modify if you've processed
% and renamed the data, but it's better to keep 'exy' as default.

[f,p] = uigetfile('E:\Ti7Al_E1_insitu_tension\Ti7Al_E1_Images\stitched_DIC','select dic files need to be rotated','multiselect','on');
if ~iscell(f)
    f = cellstr(f);
end

for iF = 1:length(f)
   load([p,f{iF}], 'exy');
   exy = -exy;
   exy_corrected = 1;
   save([p,f{iF}],'exy','exy_corrected','-append');
end