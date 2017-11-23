% chenzhe, 2017-07-14
% single pixel resolution solution for image shift.
% wrt: position of img2 with respect to img1, e.g., 'ur'(up,right),'d'(down),'l'(left) etc
% crop = [crop_row_up, crop_row_down, crop_col_left, crop_col_right]
%
% use fft to find image shift
% Ref: Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup
% and their codes
%
% chenzhe, 2017-11-20, add a Satizsky-Golay filter
%
% chenzhe, 2017-11-22, Note: in this definition, 'img2' is the 'sliding
% window'.  However, I recently found that in the mathamatical definition
% of cross correlation, F(X)G, F is the sliding window.

function [r_shift, c_shift, img2_out_ffted] = fft_register(img1_image, img2_window, wrt, crop1, crop2, filterTF)

img1 = img1_image;
img2 = img2_window;

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

% make them same size.  If img_2 is smaller, pad.  If larger, crop.
[nR,nC] = size(img1);
img = zeros(nR,nC);
img(1:min(size(img1,1),size(img2,1)),1:min(size(img1,2),size(img2,2))) = img2(1:min(size(img1,1),size(img2,1)),1:min(size(img1,2),size(img2,2)));
img2 = img;

% fft
t1 = fft2(img1);
t2 = fft2(img2);
CC = ifft2(conj(t2).*t1);
CCabs = abs(CC);
% figure;surf(CCabs,'edgecolor','none'),set(gca,'ydir','reverse');
if 1==filterTF
    filtered = CCabs;
    filtered = sgolayfilt(filtered,1,5,[],1);
    filtered = sgolayfilt(filtered,1,5,[],2);
    CCabs = CCabs-filtered;
end
[r_shift, c_shift] = find(CCabs == max(CCabs(:)),1,'first');

CCmax = CC(r_shift,c_shift)*nR*nC;

% find shift in pixel rather than index
Rs = ifftshift(-fix(nR/2):ceil(nR/2)-1);
Cs = ifftshift(-fix(nC/2):ceil(nC/2)-1);
r_shift = Rs(r_shift);
c_shift = Cs(c_shift);

% shifted image
[Cs,Rs] = meshgrid(Cs,Rs);
t2 = t2.*exp(1i*2*pi*(-r_shift*Rs/nR-c_shift*Cs/nC));
diffPhase = angle(CCmax);   % global phase diff
img2_out_ffted = abs(ifft2(t2*exp(1i*diffPhase)));

if ~exist('wrt','var')
    wrt = '';
end
if any((wrt(:)=='r')|(wrt(:)=='R'))
    if c_shift<0
        c_shift = mod(c_shift,nC);
    end
end
if any((wrt(:)=='l')|(wrt(:)=='l'))
    if c_shift>0
        c_shift = mod(c_shift,nC) - nC;
    end
end
if any((wrt(:)=='u')|(wrt(:)=='U'))
    if r_shift>0
        r_shift = mod(r_shift,nR) - nR;
    end
end
if any((wrt(:)=='d')|(wrt(:)=='D'))
    if r_shift<0
        r_shift = mod(r_shift,nR);
    end
end

r_shift = r_shift + r1a - r2a;
c_shift = c_shift + c1a - c2a;

end