% chenzhe, 2017-08-01
% similar to fft_register.  
% crop = [crop_row_up, crop_row_down, crop_col_left, crop_col_right]
%
% chenzhe, 2017-11-28, modify to use normxcorr2A()
% Can crop.  Can use filter.
%

function [yOffset, xOffSet] = normxcorr2A_register(img1_template, img2_signal, crop1, crop2, filterTF)

img1 = double(img1_template);
img2 = double(img2_signal);

r1a = 0;r1b=0;c1a=0;c1b=0;r2a=0;r2b=0;c2a=0;c2b=0;
% initial crop data
if exist('crop1','var')
    r1a = round(crop1(1));
    r1b = round(crop1(2));
    c1a = round(crop1(3));
    c1b = round(crop1(4));
    img1 = img1(1+r1a:end-r1b,1+c1a:end-c1b);
end
if exist('crop2','var')
    r2a = round(crop2(1));
    r2b = round(crop2(2));
    c2a = round(crop2(3));
    c2b = round(crop2(4));
    img2 = img2(1+r2a:end-r2b,1+c2a:end-c2b);
end


% img1 = mat2gray(img1);
% img2 = imhistmatch(mat2gray(img2),img1);
% myplot(img1);myplot(img2);

% normxcorr2
cc = normxcorr2A(img1,img2,filterTF);
[ypeak, xpeak] = find(cc == max(cc(:)));
% img1 wrt img2
yOffset = ypeak - size(img1,1);
xOffSet = xpeak - size(img1,2);

% consider initial cut
yOffset = yOffset - r1a + r2a;
xOffSet = xOffSet - c1a + c2a;

end