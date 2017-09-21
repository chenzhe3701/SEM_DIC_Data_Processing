% reduce
function gb = shrink_boundary(gb)

% gb = conv2(gb,[0 1 0; 1 1 1; 0 1 0],'same');
gb = conv2(gb,ones(3),'same');
gb(gb~=9)=0;
gb(gb~=0)=1;