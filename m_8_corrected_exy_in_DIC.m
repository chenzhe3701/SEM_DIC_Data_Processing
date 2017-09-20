
% chenzhe, 2017-06-08.
% correct the exy field of DIC data, i.e., flip the sign.
% add a note exy_corrected = 1.

[f,p] = uigetfile('D:\Marissa_test_20170430_stitched_DIC\xyuv_dic_NonIncremental','select dic files need to be rotated','multiselect','on');
if ~iscell(f)
    f = cellstr(f);
end

for iF = 1:length(f)
   load([p,f{iF}], 'exy');
   exy = -exy;
   exy_corrected = 1;
   save([p,f{iF}],'exy','exy_corrected','-append');
end