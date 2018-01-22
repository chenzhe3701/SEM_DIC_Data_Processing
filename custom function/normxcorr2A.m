% chenzhe, 2017-11-28
% similar to normxcorr2 but T can be larger than A.
% 
% reference:  
% Lewis, J. P., "Fast Normalized Cross-Correlation," Industrial Light & Magic

function c = normxcorr2A(T,A, filterTF)

% data type
T = double(T);
A = double(A);

numerator = fft_xcorr2(A,T-mean(T(:)),0);

[m, n] = size(T);
mn = m*n;
denom_T = sqrt(mn-1)*std(T(:));

local_sum_A = local_sum(A,m,n);
local_sum_A2 = local_sum(A.*A,m,n);

% raw
% local_mean_A = local_sum_A/mn;
% denom_A2 = local_sum_A2 - 2*local_sum_A*local_mean_A + (local_mean_A.^2)*mn;

% effectively, due to the relationship between sum and mean.
denom_A = sqrt( local_sum_A2 - local_sum_A.*local_sum_A/mn);

c = numerator./denom_T./denom_A;
c(0 == denom_A) = 0;

if (exist('filterTF','var') && (1==filterTF))
    filtered = c;
    filtered = sgolayfilt(filtered,1,5,[],1);
    filtered = sgolayfilt(filtered,1,5,[],2);
    c = c-filtered;
end

end


% imaging we have a 2d-cumsum [AB;CD], then D = ABCD-AC-AB+A
% [AB;CD] is 'fen kuai ju zhen'
function s = local_sum(A,m,n)
B = A;
B = padarray(B,[m,n],'pre');
B = padarray(B,[m-1,n-1],'post');
B = cumsum(B,1);
B = cumsum(B,2);
s = B(1+m:end,1+n:end) - B(1+m:end, 1:end-n) - B(1:end-m, 1+n:end) + B(1:end-m, 1:end-n);
end