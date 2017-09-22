% [tform, A] = make_avg_affine(cpFrom, cpTo)
% make an averaged version of affine transformation
% [cpTo,1] = A * [cpFrom';1]
% tform = maketform('affine',A' );  ---> note it's A prime.
% chenzhe, 2017-07-10
%
% chenzhe, 2017-09-21
% A is my usually way of thinking, points on right
% T = A', maketform('affine', T)
% Here according to definition, cpTo(n x (2+1)) = cpFrom(n x (2+1)) x A' 
% Q = PA', or Q' = AP'
% A' = [a b 0; c d 0; e f 1].
% A' = pinv(P) * Q, or A = Q' * pinv(P')
% P is col-full-rank, so inv(P) should use left, inv(P') use right.

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

% % equivalently, Q=PT, so T=inv(P) * Q
% T = inv(P'*P)*P'*Q;
% T(1:end-1,3)=0;
% T(end,3) = 1;
% tform = maketform('affine',T );