% chenzhe, 2017-12-07
% use the definition
%
% T = [h1,h4,h7;
%      h2,h5,h8;
%      h3,h6,h9];
% 
% for 'affine', h7=h8=0, h9=1
% [x_hat, y_hat, 1] = [x,y,1] * T;
%  
% for 'projective', h9=1
% x_hat = [x,y,1]*T(:,1) / [x,y,1]*T(:,3)
% y_hat = [x,y,1]*T(:,2) / [x,y,1]*T(:,3)
%
% the 'affine' can actually use the generalized form in 'projective'  
%
% this code solves a system of equations, by writing into a vector, then
% 0 = [x1,y1,1, 0,0,0, -x1_hat*x1,-x1_hat*y1,-x1_hat*1] * h[9x1]
% 0 = [0,0,0, x1,y1,1, -y1_hat*x1,-y1_hat*y1,-y1_hat*1] * h[9x1]
% 0 = [x2,y2,1, 0,0,0, -x2_hat*x2,-x2_hat*y2,-x2_hat*1] * h[9x1]
% 0 = [0,0,0, x2,y2,1, -y2_hat*x2,-y2_hat*y2,-y2_hat*1] * h[9x1]
% ...
%
% then, pick 6 or 8 unknowns to solve, depending on type.

function [tform, T] = make_average_transform(type, cpFrom, cpTo)

x_hat = cpTo(:,1);
y_hat = cpTo(:,2);

xy1 = [cpFrom,ones(size(cpFrom,1),1)];

k = [];
for ii=1:size(x_hat,1)
   k = [k; 
       kron([1 0 -x_hat(ii)], xy1(ii,:));
       kron([0 1 -y_hat(ii)], xy1(ii,:));
       ];
end

if strcmpi(type, 'affine')
   A = k(:,1:6);
   B = -k(:,9);
   T = A\B;
   % T = pinv(A'*A)*A'*B;
   T = reshape([T(:);0;0;1],3,3);
   tform = maketform('affine',T);
elseif strcmpi(type, 'projective')
   A = k(:,1:8);
   B = -k(:,9);
   T = A\B;
   % T = pinv(A'*A)*A'*B;
   T = reshape([T(:);1],3,3);
   tform = maketform('projective',T);
end


