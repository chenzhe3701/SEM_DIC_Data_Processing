% Chenzhe, 2017-05-12
% It's weird, but sometimes you know you are not looking at the sample, you
% still get signal for ETD and DIC gives you measurement.  So I'd like to
% zeros these sigma values.
%
% append a variable called 'mask'

%% select files

[f,p] = uigetfile('D:\Marissa_test_20170430_FlashDriveData\ts5Al_02_all images_sorted FOV _non_incremental_ref_e1\','Select .mat fileS of the same FOV','multiselect','on');
if ~iscell(f)
    f = cellstr(f);
end

%% creat a maks

%  load([p,'SigmaMask'],'mask');

data = load([p,f{end}],'sigma');
imagesc(data.sigma); title('Draw polygon over good region. Double click to confirm. Then X to quit.');
H = impoly;
H.wait;
disp('mask created');
mask = H.createMask;
imagesc(mask);title('This is the mask of sigma not zero');
mask = double(mask);

%% If you are satisfield with the mask
% temporarily just for these few fields.

% This append a variable called mask.
AppendMask = 1;
if(AppendMask==1)
    for ii = 1:length(f)
        save([p,f{ii}],'mask','-append');
    end
end


% % This correct the data. Not favored.
% for ii = 1:length(f)
%     clear u v x y exx exy eyy e1 e2 U V X Y
%     copyfile([p,f{ii}],[p,'old_',f{ii}]);
%     load([p,f{ii}]);
%     sigma(~mask) = -1;
%     u(~mask) = nan;
%     v(~mask) = nan;
%     exx(~mask) = nan;
%     exy(~mask) = nan;
%     eyy(~mask) = nan;
%     save([p,f{ii}],'exx','exy','eyy','sigma','u','v','x','y','-append');
% end

save([p,'SigmaMask'],'mask');


