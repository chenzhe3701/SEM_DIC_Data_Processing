
% chenzhe, 2017-11-19
% using matlab fft to calculate xcorrelation.
%
% by doing this padding, I believe this gives the same result as xcorr2,
% but I force it to used the faster method.


function c = fft_xcorr2(A,B)
[r1,c1] = size(A);
[r2,c2] = size(B);

AA = zeros(r1+r2-1, c1+c2-1);
AA(r2:end,c2:end) = A;

BB = zeros(r1+r2-1, c1+c2-1);
BB(1:r2, 1:c2)=B;

c = ifft2(fft2(AA).*conj(fft2(BB)));
end