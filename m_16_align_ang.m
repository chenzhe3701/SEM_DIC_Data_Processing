% try to align ang/ebsd data.  This is more of a manual mode.
% 
% based on some fft-based image align codes and thoughts
% 
% chenzhe, 2017-07-14
% chenzhe, 2017-11-27, update to use filter

f1 = 'd:\Ti7Al_N3_post1_left.ang';
f2 = 'd:\Ti7Al_N3_post1_right.ang';
[p1,p,p2,x,y,IQ1] = read_ang(f1);
[P1,P,P2,~,~,IQ2] = read_ang(f2);

% % cut img_1 from upper_left, and img_2 from lower_right
% rCut = 0;
% cCut = 0;
% % additional cut of img_2 from upper_left, and img_1 also from upper_left
% rCut2 = 0;
% cCut2 = 0;
% map1 = zeros(size(IQ1)-[rCut+rCut2,cCut+cCut2]);
% map2 = zeros(size(IQ1)-[rCut+rCut2,cCut+cCut2]);
% for iR = 1:size(map1,1)
%    for iC = 1:size(map1,2)
%       map1(iR,iC) = calculate_misorientation_euler_d([0 0 0], [p1(iR+rCut+rCut2,iC+cCut+cCut2),p(iR+rCut+rCut2,iC+cCut+cCut2),p2(iR+rCut+rCut2,iC+cCut+cCut2)]/pi*180,'hcp');
%       map2(iR,iC) = calculate_misorientation_euler_d([0 0 0], [P1(iR+rCut2,iC+cCut2),P(iR+rCut2,iC+cCut2),P2(iR+rCut2,iC+cCut2)]/pi*180,'hcp');
%    end
%    disp(iR);
% end

% This part is very slow.
map1 = cell2mat(arrayfun(@(x,y,z) calculate_misorientation_euler_d([0,0,0],[x,y,z]/pi*180,'hcp'), p1, p, p2,'uniformoutput',false));
map2 = cell2mat(arrayfun(@(x,y,z) calculate_misorientation_euler_d([0,0,0],[x,y,z]/pi*180,'hcp'), P1, P, P2,'uniformoutput',false));


% [r_shift,c_shift] = fft_register(map1,map2,'r',[0,0, 600,0],[0,0, 0,600],0) % without filter, does not work well
[r_shift,c_shift] = fft_register(map1,map2,'r',[0,0, 0,0],[0,0, 0,0],1)
figure;
subplot(1,3,1);
imagesc(map1);
subplot(1,3,2);
imagesc(map2);
subplot(1,3,3);
imagesc(circshift(circshift(map2,r_shift,1),c_shift,2));

% this is a double check
[r_shift_2,c_shift_2] = fft_register(IQ1,IQ2,'r',[0,0, 0,0],[0,0, 0,0],1)
figure;
subplot(1,3,1);
imagesc(IQ1);
subplot(1,3,2);
imagesc(IQ2);
subplot(1,3,3);
imagesc(circshift(circshift(IQ2,r_shift_2,1),c_shift_2,2));
