% function pooled = my_pool(img, poolSize, method)
% 
% chenzhe 2017-04-10
% a function for local pooling
% method: 1=mean, 2=max

function pooled = my_pool(img, poolSize, method)

[nR,nC] = size(img);
nR_pool = floor(nR/poolSize);
nC_pool = floor(nC/poolSize);
pooled = zeros(nR_pool, nC_pool);
for iR = 1:nR_pool
    for iC = 1:nC_pool
        patch = img(poolSize*(iR-1)+1:poolSize*iR,poolSize*(iC-1)+1:poolSize*iC);
        if method == 1
            pooled(iR,iC) = mean(patch(:));
        elseif method == 2
            pooled(iR,iC) = max(patch(:));
        end
    end
end

end