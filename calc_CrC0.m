function [C_0,C_r,mu_i] = calc_CrC0(curr_frame,next_frame,mode,tou_c,gamma)

global offset_x;
global offset_y;
global f_count_f;

% mode 0 : calc C_0,C_r only
% mode 1 : calc mu_i additionally (only used when frame t,t+1)

curr_frame_gray = rgb2gray(curr_frame);
next_frame_gray = rgb2gray(next_frame);

[f_height,f_width] = size(curr_frame_gray);
%{
  autodetected harris position is better than using fix # of features
%}
%harris_curr = corner(curr_frame_gray,500); 
%harris_next = corner(next_frame_gray,500);
harris_curr = detectHarrisFeatures(curr_frame_gray);
harris_next = detectHarrisFeatures(next_frame_gray);
%imshow(curr_frame_gray), hold on,plot(harris_curr(:,1), harris_curr(:,2), 'r*');
%figure,imshow(next_frame_gray), hold on,plot(harris_next(:,1), harris_next(:,2), 'r*');

[feature_curr,valid_curr] = extractFeatures(curr_frame_gray,harris_curr,'Method','FREAK');
[feature_next,valid_next] = extractFeatures(next_frame_gray,harris_next,'Method','FREAK');
% imshow(curr_frame_gray),title('Base image'),hold on,plot(valid_curr);
% figure; imshow(next_frame_gray);title('Transformed image');hold on, plot(valid_next);

index_pairs = matchFeatures(feature_curr,feature_next,'Unique',true); % exhaustive match, 1-1 match

match_curr  = valid_curr(index_pairs(:,1),:);
match_next = valid_next(index_pairs(:,2),:);

%{
[tform,inlier_1,inlier_2] = estimateGeometricTransform(match_1,match_2,transformType)
 tform maps the inliers in match_1 to the inliers in match_2.
 tform.T having the convention below
 [x y 1] = [u v 1] * T
    where T has the form:
        [a b 0;
        c d 0;
        e f 1];
%}
[tform,inlier_next,inlier_curr,status] = estimateGeometricTransform(match_next,match_curr,'affine', ...
                                         'Confidence',99,'MaxDistance',sqrt(f_height^2+f_width^2));
if status ~= 0 
    % can't find transform, in this case, believe that Cr>tou_c, so set C0
    % to gamma, and Cr to inf
    C_0 = gamma;
    C_r = 100000;
    mu_i = NaN;
    return
end

inlier_next_coor = inlier_next.Location;
inlier_curr_coor = inlier_curr.Location;
T = tform.T;

warping_result = padarray(inlier_next_coor,[0 1],1,'post') * T;
warping_result = (warping_result(:,1:2) - inlier_curr_coor).^2;
warping_result = sqrt(warping_result(:,1)+warping_result(:,2));
C_r = sum(warping_result) / numel(warping_result);

[h,w] = size(curr_frame_gray);
img_center = [w/2 h/2 1];
warping_result = img_center * T;
warping_result = (warping_result - img_center).^2;
if C_r <= tou_c
    C_0 = sqrt(warping_result(:,1)+warping_result(:,2));
else
    %{
    below case do exist... possibly due to motion blur even with
    consecutive time
    if mode == 1
        error('consecutive but with Cr>tou_c \n')
    end
    %}
    C_0 = gamma;
end

if mode==1
    img_boundary = [1 h/2 1; w/2 1 1;
                    1  1  1;  w  h 1;];
    warping_result = img_boundary * T;
    warping_result = (warping_result - img_boundary).^2;
    mu_i = sum(sqrt(warping_result(:,1)+warping_result(:,2)))/4;
    offset_x = [offset_x T(3,1)];
    offset_y = [offset_y T(3,2)];
    f_count_f = f_count_f + 1;
elseif mode==0
    mu_i = 0;
else 
    fprintf('error mode ! \n');
end

%figure; showMatchedFeatures(curr_frame_gray,next_frame_gray,match_curr,match_next);
%title('Matched SURF points,including outliers');
%figure; showMatchedFeatures(curr_frame_gray,next_frame_gray,inlier_curr,inlier_next);
%title('Matched inlier points');
