% Normalise the points to the centroid
function [newpts, T] = NormalisePoints(pts)
    % pts is a 3 row vector
    if size(pts,1) ~= 3
        error('pts must be 3xN');
    end
    
    % Find the indices of the points that are not at infinity
    % eps is the smallest number in matlab, which means 3rd row can't be
    % zero
    finiteind = find(abs(pts(3,:)) > eps); 
    
    if length(finiteind) ~= size(pts,2)
        disp('Some points are at infinity');
    end
    
    % For the finite points ensure homogeneous coords have scale of 1
    pts(1,finiteind) = pts(1,finiteind)./pts(3,finiteind);
    pts(2,finiteind) = pts(2,finiteind)./pts(3,finiteind);
    pts(3,finiteind) = 1;
    
    % mean calculate along column, so first need to convert pts to col
    % vector
    c = mean(pts(1:2,finiteind)')';  
    newp(1,finiteind) = pts(1,finiteind)-c(1); 
    newp(2,finiteind) = pts(2,finiteind)-c(2);
    
    dist = sqrt(newp(1,finiteind).^2 + newp(2,finiteind).^2);
    meandist = mean(dist(:)); 
    
    scale = sqrt(2)/meandist;
    
    T = [scale   0   -scale*c(1)
         0     scale -scale*c(2)
         0       0      1      ];
    
    newpts = T*pts;