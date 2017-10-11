function warpIm = myWarp(minx,maxx,miny,maxy,im,warpIm,invH,gap)

    [h,w,~] = size(warpIm);
    
    if gap > 0
        minx = floor(minx);
        miny = floor(miny);
        
        % added to prevent largest offset+gap < 0
        if minx + gap <= 0 
            minx = -gap+1;
        end
        if miny + gap <= 0
            miny = -gap+1;
        end
        %minx = max(minx,1);
        %miny = max(miny,1);
        
    else
        minx = max(floor(minx),1);
        miny = max(floor(miny),1);
        
        % possible minx or miny > w or h (happens when track is long,
        % warping too large)
        minx = min(minx,w);
        miny = min(miny,h);
    end
    
    
    
    maxx = min(ceil(maxx),w);
    maxy = min(ceil(maxy),h);

    %szIm = [maxy-miny,maxx-minx];
    szIm = [maxy-miny+1,maxx-minx+1];
    %[x,y] = meshgrid(minx:maxx-1,miny:maxy-1);
    [x,y] = meshgrid(minx:maxx,miny:maxy);
    pix   = [x(:)'; y(:)']; %  convert x,y to homogeneous coordinate representation
    hPixels = [ pix; ones(1,prod(szIm))];
    hScene  = invH*hPixels;
    xprime=(hScene(1,:)./(hScene(3,:)))';
    yprime=(hScene(2,:)./(hScene(3,:)))';
    
    xprime = reshape(xprime, szIm); % warped coordinate 
    yprime = reshape(yprime, szIm);
    
    mask_x = (xprime>=1) .* (xprime<=w);
    mask_y = (yprime>=1) .* (yprime<=h);
    mask = logical(mask_x .* mask_y);
    
    result(:,:,1) = interp2(im(:,:,1),xprime,yprime,'cubic');
    result(:,:,2) = interp2(im(:,:,2),xprime,yprime,'cubic');
    result(:,:,3) = interp2(im(:,:,3),xprime,yprime,'cubic');
    
%   %use spline interpolation for better quality with lower speed
%   result(:,:,1) = interp2(im(:,:,1),xprime,yprime,'spline');
%   result(:,:,2) = interp2(im(:,:,2),xprime,yprime,'spline');
%   result(:,:,3) = interp2(im(:,:,3),xprime,yprime,'spline');
    
    result_tmp = result(:,:,1);
    tmp = warpIm(:,:,1);
    tmp(mask) = result_tmp(mask);
    warpIm(:,:,1) = tmp;
    result_tmp = result(:,:,2);
    tmp = warpIm(:,:,2);
    tmp(mask) = result_tmp(mask);
    warpIm(:,:,2) = tmp;
    result_tmp = result(:,:,3);
    tmp = warpIm(:,:,3);
    tmp(mask) = result_tmp(mask);
    warpIm(:,:,3) = tmp;
    %warpIm(miny+gap:maxy-1+gap,minx+gap:maxx-1+gap,1) = result(:,:,1);
    %warpIm(miny+gap:maxy-1+gap,minx+gap:maxx-1+gap,2) = result(:,:,2);
    %warpIm(miny+gap:maxy-1+gap,minx+gap:maxx-1+gap,3) = result(:,:,3);

end