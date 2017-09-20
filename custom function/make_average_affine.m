% [tform, A] = make_avg_affine(cpFrom, cpTo)
% make an averaged version of affine transformation
% [cpTo,1] = A * [cpFrom';1]
% tform = maketform('affine',A');
% chenzhe, 2017-07-10

function [tform, A] = make_average_affine(cpFrom, cpTo)
if (size(cpFrom,2) ~= 2)||(size(cpTo,2)~=2)
    error('input points make n by 2 matrices');
end
% Q' = A*P';   note, a more mathamatical way should be Q = PA, A on right
Q = [cpTo, ones(size(cpFrom,1),1)];
P = [cpFrom, ones(size(cpTo,1),1)];
A = Q'*P*inv(P'*P);   % averaged affine
A(3,1:end-1)=0;
A(3,end) = 1;
tform = maketform('affine',A');