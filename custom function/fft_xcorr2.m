
% chenzhe, 2017-11-19
% using matlab fft to calculate xcorrelation.
%
% by doing this padding, I believe this gives the same result as xcorr2,
% but I force it to used the faster method.
%
% 'F' is the filter.
% Math definition, corresponding to F(xcorr)G.
% Physical interpretation, F sldes on G.  So, upperleft element in 'c'
% means F is upper-left with respect to G.


function c = fft_xcorr2(G_img_signal, F_filter_window, filterTF)
[r1,c1] = size(G_img_signal);
[r2,c2] = size(F_filter_window);

explicit = 0;
if explicit
    % explicit padding
    
    G = zeros(r1+r2-1, c1+c2-1);
    G(r2:end,c2:end) = G_img_signal;
    
    F = zeros(r1+r2-1, c1+c2-1);
    F(1:r2, 1:c2) = F_filter_window;
    
    c = ifft2(conj(fft2(F)).*fft2(G));
else
    % effectively, padding
    c = circshift( ifft2( conj(fft2(F_filter_window, r1+r2-1, c1+c2-1)) .* fft2(G_img_signal, r1+r2-1, c1+c2-1) ), [(r2-1), (c2-1)] );
end

% if use filter
if (exist('filterTF','var') && (1==filterTF))
    filtered = c;
    filtered = sgolayfilt(filtered,1,5,[],1);
    filtered = sgolayfilt(filtered,1,5,[],2);
    c = c-filtered;
end

end