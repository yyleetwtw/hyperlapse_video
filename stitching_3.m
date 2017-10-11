function final_img =stitching_3(frame_0,frame_1,frame_2,T_0to1,T_2to1,B_1opt)

%{
 first decode if frame_0, frame_2 is usable, frame_1 is current frame,
 this function stitching f0,f2 to f1 with an eye to reducing temporal
 copping area
%}

f0_usable = ~isempty(find(frame_0,1));
f2_usable = ~isempty(find(frame_2,1));

[fheight,fwidth,no_use] = size(frame_1);

src = [1 fwidth 1 fwidth;
       1 1  fheight fheight;
       1 1 1 1;];
%B_1opt = P_queue(:,:,1)*inv(buf_queue(1).F);
B_1opt = B_1opt ./ B_1opt(3,3);
dst_1 = B_1opt * src;
if f0_usable
    H_0to1opt = B_1opt*T_0to1;
    H_0to1opt = H_0to1opt ./ H_0to1opt(3,3);
else
    H_0to1opt = eye(3);
end
dst_0to1 = H_0to1opt * src;
if f2_usable
    H_2to1opt = B_1opt*T_2to1;
    H_2to1opt = H_2to1opt ./ H_2to1opt(3,3);
else
    H_2to1opt = eye(3);
end
dst_2to1 = H_2to1opt * src;

% Find the minimum and maximum output limits
xMin = min([1 dst_0to1(1,:) dst_2to1(1,:) dst_1(1,:)]);
xMax = max([fwidth dst_0to1(1,:) dst_2to1(1,:) dst_1(1,:)]);

yMin = min([1 dst_0to1(2,:) dst_2to1(2,:) dst_1(1,:)]);
yMax = max([fheight 1 dst_0to1(2,:) dst_2to1(2,:) dst_1(1,:)]);

% Width and height of panorama.
width  = round(xMax - xMin);
height = round(yMax - yMin);
if width < fwidth
    width = fwidth;
end
if height < fheight
    height = fheight;
end

% Initialize the "empty" panorama.
panorama = zeros([height width 3], 'like', frame_1);

blender = vision.AlphaBlender('Operation', 'Binary mask', ...
    'MaskSource', 'Input port');

% Create a 2-D spatial reference object defining the size of the panorama.
xLimits = [xMin xMax];
yLimits = [yMin yMax];
panoramaView = imref2d([height width], xLimits, yLimits);

if f0_usable
    % Transform frame_0 into the panorama.
    tmp = projective2d(eye(3));
    %H_0to1opt(3,1:2) = 0;
    tmp.T = H_0to1opt';
    warpedImage = imwarp(frame_0, tmp, 'OutputView', panoramaView);
    panorama = step(blender, panorama, warpedImage, warpedImage(:,:,1));
end
if f2_usable
    % Transform frame_2 into the panorama.
    %H_2to1opt(3,1:2) = 0;
    tmp.T = H_2to1opt';
    warpedImage = imwarp(frame_2, tmp, 'OutputView', panoramaView);
    panorama = step(blender, panorama, warpedImage, warpedImage(:,:,1));
end

B_1opt(3,1:2) = 0;
tmp.T = B_1opt';
warpedImage = imwarp(frame_1, tmp, 'OutputView', panoramaView);
panorama = step(blender, panorama, warpedImage, warpedImage(:,:,1));

mask_x = (xMin<=0);
mask_y = (yMin<=0);
orig_pos = [abs(xMin)*mask_x+1 abs(xMin)*mask_x+fwidth; abs(yMin)*mask_y+1 abs(yMin)*mask_y+fheight];
orig_pos = floor(orig_pos);
orig_pos(1,1) = (1-mask_x)*1 + mask_x*orig_pos(1,1);
orig_pos(1,2) = (1-mask_x)*fwidth + mask_x*orig_pos(1,2);
orig_pos(2,1) = (1-mask_y)*1 + mask_y*orig_pos(2,1);
orig_pos(2,2) = (1-mask_y)*fheight + mask_y*orig_pos(2,2);

final_img = zeros(fheight,fwidth,3);
final_img(:,:,1) = panorama(orig_pos(2,1):orig_pos(2,2), orig_pos(1,1):orig_pos(1,2),1);
final_img(:,:,2) = panorama(orig_pos(2,1):orig_pos(2,2), orig_pos(1,1):orig_pos(1,2),2);
final_img(:,:,3) = panorama(orig_pos(2,1):orig_pos(2,2), orig_pos(1,1):orig_pos(1,2),3);



%{
buildingDir = fullfile(toolboxdir('vision'), 'visiondata', 'building');
buildingScene = imageSet(buildingDir);
montage(buildingScene.ImageLocation)
I = read(buildingScene, 1);
grayImage = rgb2gray(I);
points = detectSURFFeatures(grayImage);
[features, points] = extractFeatures(grayImage, points);
tforms(buildingScene.Count) = projective2d(eye(3));
for n = 2:buildingScene.Count
    
    % Store points and features for I(n-1).
    pointsPrevious = points;
    featuresPrevious = features;
    
    % Read I(n).
    I = read(buildingScene, n);
    
    % Detect and extract SURF features for I(n).
    grayImage = rgb2gray(I);
    points = detectSURFFeatures(grayImage);
    [features, points] = extractFeatures(grayImage, points);
    
    % Find correspondences between I(n) and I(n-1).
    indexPairs = matchFeatures(features, featuresPrevious, 'Unique', true);
    
    matchedPoints = points(indexPairs(:,1), :);
    matchedPointsPrev = pointsPrevious(indexPairs(:,2), :);
    
    % Estimate the transformation between I(n) and I(n-1).
    tforms(n) = estimateGeometricTransform(matchedPoints, matchedPointsPrev,...
        'projective', 'Confidence', 99.9, 'MaxNumTrials', 2000);
    
    % Compute T(1) * ... * T(n-1) * T(n)
    tforms(n).T = tforms(n-1).T * tforms(n).T;
end
imageSize = size(I);  % all the images are the same size
for i = 1:numel(tforms)
    [xlim(i,:), ylim(i,:)] = outputLimits(tforms(i), [1 imageSize(2)], [1 imageSize(1)]);
end
avgXLim = mean(xlim, 2);
[~, idx] = sort(avgXLim);
centerIdx = floor((numel(tforms)+1)/2);
centerImageIdx = idx(centerIdx);
Tinv = invert(tforms(centerImageIdx));
for i = 1:numel(tforms)
    tforms(i).T = Tinv.T * tforms(i).T;
end
for i = 1:numel(tforms)
    [xlim(i,:), ylim(i,:)] = outputLimits(tforms(i), [1 imageSize(2)], [1 imageSize(1)]);
end
xMin = min([1; xlim(:)]);
xMax = max([imageSize(2); xlim(:)]);
yMin = min([1; ylim(:)]);
yMax = max([imageSize(1); ylim(:)]);
width  = round(xMax - xMin);
height = round(yMax - yMin);
panorama = zeros([height width 3], 'like', I);
blender = vision.AlphaBlender('Operation', 'Binary mask', ...
    'MaskSource', 'Input port');
xLimits = [xMin xMax];
yLimits = [yMin yMax];
panoramaView = imref2d([height width], xLimits, yLimits);
for i = 1:buildingScene.Count   
    I = read(buildingScene, i);
    % Transform I into the panorama.
    warpedImage = imwarp(I, tforms(i), 'OutputView', panoramaView);
    % Overlay the warpedImage onto the panorama.
    panorama = step(blender, panorama, warpedImage, warpedImage(:,:,1));
end
figure,imshow(panorama)
%}