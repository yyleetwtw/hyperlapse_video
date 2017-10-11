clear all,clc,close all

% modified external library
addpath('MeshModel\mesh');
addpath('MeshModel\RANSAC');

% ----- mode
WRITE_WARP = 1;
DRAW_PATH = 0;

% ----- get input from user
prompt = 'What is the selected video file name? ';
selected_frame_filename = input(prompt,'s')
prompt = 'What is the hyperlapse video file name? ';
Hyperlapse_name = input(prompt,'s')
%prompt = 'What is smooth weight (0.5~1)? ';
%smooth_weight = input(prompt)
smooth_weight = 0.75;

% ----- parameters
speedup_rate = 10;
video_seq=VideoReader(selected_frame_filename);
Hyperlapse = VideoWriter(Hyperlapse_name,'MPEG-4');
fwidth = video_seq.Width;
fheight = video_seq.Height;
num_frame = floor(video_seq.FrameRate*video_seq.Duration); % count from 1
transform_thres = sqrt(fwidth^2+fheight^2)*5/8;
cummulate_thres = 500;
gau_std = sqrt(fwidth^2+fheight^2)/8;
w_tr_window = 30/speedup_rate;

track_cnt = 1;
track_cummulate = 0;
%F_queue = []; % used to store struct array of mat
queue_cnt = 1;
all_trans_dist = [];
fc_seq = [];

%% compute mesh model
fc = 1;
fprintf('Computing mesh model start\n');
while fc < num_frame

    if fc <= 1
        frame_curr = readFrame(video_seq);
        fc = fc+1;
        continue;
    end
    
    frame_prev = frame_curr;
    if hasFrame(video_seq)
        frame_curr = readFrame(video_seq);
        [fc]
    else
        fprintf('No frame can be read');
        break
    end
    
    [src_inlier,target_inlier,state]=SURF(frame_prev,frame_curr);  % ret pair (src, target)
    if state == 1
        fprintf('lost track of feature, reinitialize a track \n');
        if track_cummulate ~= 0 % case that not 2 consecutive reinitialize
            track_cnt = track_cnt + 1;
        end
        track_cummulate = 0;
        fc = fc + 1;
        continue;
    else
        fprintf('[DONE]\n');
    end
    
    if length(src_inlier) < 2
        fprintf('not enough matched features, reinitialize a track \n');
        if track_cummulate ~= 0 % case that not 2 consecutive reinitialize
            track_cnt = track_cnt + 1;
        end
        track_cummulate = 0;
        fc = fc + 1;
        continue;
    end
    
    if track_cummulate > cummulate_thres  % limit cummulate # to prevent undesired result
        fprintf('cummulate limite arrived, reinitialize a track \n');
        track_cummulate = 0;
        track_cnt = track_cnt + 1;
        continue;
    end

    quadWidth = fwidth/(1);
    quadHeight = fheight/(1);

    num_feature = numel(src_inlier)/2;
    lamda = 0.3; %mesh more rigid if larger value. [0.2~5], corresponding to alpha in paper
    first_alpha = 0;
    F = zeros(3,3);
    while lamda <= 3
        asap = AsSimilarAsPossibleWarping(fheight,fwidth,quadWidth,quadHeight,lamda);
        asap.SetControlPts(src_inlier,target_inlier);%set matched features
        asap.Solve();            %solve Ax=b for as similar as possible
        homos = asap.CalcHomos();% calc local hommograph transform
        for j=1:3
            for i=1:3
                F(j,i) = homos(1,1,j,i); % the (1,1) only for 1x1 cell grid
            end
        end

        src_transformed = F * [src_inlier'; ones(1,num_feature)];
        for col = 1:num_feature
            src_transformed(:,col) = src_transformed(:,col)./src_transformed(3,col);
        end
        src_transformed = src_transformed';
        src_transformed = src_transformed(:,1:2);
        transform_dist = sum(sum((target_inlier - src_transformed).^2))/num_feature;
        if first_alpha == 0
            curr_min_dist = transform_dist;
            F_best = F;
            first_alpha = 1;
            best_alpha = lamda;
        else 
            if curr_min_dist > transform_dist
                curr_min_dist = transform_dist;
                F_best = F;
                best_alpha = lamda;
            end
        end
        lamda = lamda + 0.3;
    end
    
    all_trans_dist = [all_trans_dist curr_min_dist];
    
    if curr_min_dist > transform_thres
        fprintf('too large transform distance, reinitialize a track \n');
        if track_cummulate ~= 0 % case that not 2 consecutive reinitialize
            track_cnt = track_cnt + 1;
        end
        track_cummulate = 0;
        fc = fc + 1;
        continue;
    end
    
    track_cummulate = track_cummulate+1;
    tmp.F = F;
    tmp.track = track_cnt;
    tmp.track_cum = track_cummulate;
    tmp.fc = fc;
    F_queue(queue_cnt) = tmp;
    queue_cnt = queue_cnt + 1;
        
    fc = fc+1;
end

[~,length_queue] = size(F_queue);
for i=1:length_queue
    fc_seq = [fc_seq F_queue(i).fc];
end

trans_dist = [];
fc_seq = [];
fprintf('Warping start\n');
[~,length_queue] = size(F_queue);

for i=1:length_queue
    F= F_queue(i).F;
    trans_dist =[trans_dist sqrt(sum(F(:,3).^2)-1)];
    fc_seq = [fc_seq F_queue(i).fc];
end

total_track = F_queue(end).track;
start_end_of_track = zeros(total_track,4); % (start_fc, end_fc, start_queue_row, end_queue_row)

for curr_track=1:total_track
    for track_row_num = 1:length_queue
        if F_queue(track_row_num).track == curr_track
            start_end_of_track(curr_track,1) = F_queue(track_row_num).fc - F_queue(track_row_num).track_cum;
            start_end_of_track(curr_track,2) = F_queue(track_row_num).fc;
        end
    end
end
for curr_track=1:total_track
    for track_row_num = 1:length_queue
        if F_queue(track_row_num).fc == start_end_of_track(curr_track,1)+1
            start_end_of_track(curr_track,3) = track_row_num;
        end
        if F_queue(track_row_num).fc == start_end_of_track(curr_track,2)
            start_end_of_track(curr_track,4) = track_row_num;
        end
    end
end

buf_queue = F_queue;

% warp using smoothed path
% first pass : generate C(t)
prev_itr_end_pos = 1;
for curr_track=1:total_track
    start_fc = start_end_of_track(curr_track,1);
    end_fc = start_end_of_track(curr_track,2);
    track_seq = [];
    for track_row_num = prev_itr_end_pos:length_queue
        if F_queue(track_row_num).fc - start_fc == 1 % first F of track, maintain the same
            % buf_queue.F corresponding to C in paper
            %buf_queue(track_num).F = buf_queue(track_num).F;
        else
            buf_queue(track_row_num).F = buf_queue(track_row_num).F * buf_queue(track_row_num-1).F;
        end
        if F_queue(track_row_num).fc == end_fc
            prev_itr_end_pos = track_row_num +1;
            break
        end
    end
end


Hyperlapse.FrameRate = 30;
open(Hyperlapse);
video_seq=VideoReader(selected_frame_filename);
fc = 0;
% grid vertex coordinate only for 1x1 cell
src = [1 fwidth 1 fwidth;
    1 1  fheight fheight;
    1 1 1 1;];
warp_img = [];
gap = 0;
P_queue = zeros(3,3,length_queue);
src_frame = zeros(video_seq.Height,video_seq.Width,3);      % curr frame
src_frame_prev = zeros(video_seq.Height,video_seq.Width,3); % prev frame
src_frame_next = zeros(video_seq.Height,video_seq.Width,3); % next frame

for curr_track=1:total_track
    % idx used to specify member position of certain track
    start_queue_row = start_end_of_track(curr_track,3);
    end_queue_row = start_end_of_track(curr_track,4);
    
    %% generate w_tr for every track
    
    w_tr = zeros(end_queue_row-start_queue_row+1 , end_queue_row-start_queue_row); % row:curr_t, col:r~=t & |t-r|<=30
    for curr_track_mem_row = 1:end_queue_row-start_queue_row+1
        w_tr_col = 1;
        C_t = buf_queue(curr_track_mem_row+start_queue_row-1).F;
        for pair_mem = 1:end_queue_row-start_queue_row+1
            if pair_mem == curr_track_mem_row
                continue
            end
            C_r = buf_queue(pair_mem+start_queue_row-1).F;
            camera_motion = abs(C_t(1,3)-C_r(1,3))+abs(C_t(2,3)-C_r(2,3));
            time_dist = buf_queue(curr_track_mem_row+start_queue_row-1).fc - buf_queue(pair_mem+start_queue_row-1).fc;
            if abs(time_dist) > w_tr_window
                w_tr(curr_track_mem_row, w_tr_col) = 0;
            else
                %w_tr(curr_track_mem, w_tr_col) = normpdf(abs(time_dist),0,10) * normpdf(camera_motion,0,gau_std);
                w_tr(curr_track_mem_row, w_tr_col) = normpdf(abs(time_dist),0,10);
            end
            w_tr_col = w_tr_col+1;
        end
    end
    [w_tr_row,~]=size(w_tr);
    for m=1:w_tr_row
        row_sum = sum(w_tr(m,:));
        if row_sum > 0
            w_tr(m,:)=w_tr(m,:)/row_sum*smooth_weight;
        end
    end
    
    %% optimize each track and obtain optimal path P
    lambda_t = 5;
    P = zeros(3,3,end_queue_row-start_queue_row+1); % (homo,homo,t)
    % initialize P(t) as C(t)
    P_prev = zeros(3,3,end_queue_row-start_queue_row+1);
    for curr_opt_mem = 1 : end_queue_row-start_queue_row+1
        P_prev(:,:,curr_opt_mem) = buf_queue(curr_opt_mem+start_queue_row-1).F;
    end
    % iterative optimization process
    opt_itr = 1;
    while opt_itr < 100
        for curr_opt_mem = 1 : end_queue_row-start_queue_row+1
            w_tr_col = 1;
            gamma = 1+2*lambda_t*sum(w_tr(curr_opt_mem,:));
            tmp = zeros(3,3);
            for pair_mem = 1:end_queue_row-start_queue_row+1
                if pair_mem == curr_opt_mem
                    continue
                end
                w = w_tr(curr_opt_mem, w_tr_col);
                tmp = tmp + 2*lambda_t*w/gamma * P_prev(:,:,pair_mem);
                w_tr_col = w_tr_col+1;
            end
            P(:,:,curr_opt_mem) = 1/gamma*buf_queue(curr_opt_mem+start_queue_row-1).F + tmp;
        end
        P_prev = P;
        opt_itr = opt_itr + 1;
    end
    P_queue(:,:,start_queue_row:end_queue_row) = P;
    
    %% below warping process
    % matrices used to visualize opt path & orig path
    camera_x = [];
    camera_y = [];
    opt_x = [];
    opt_y = [];
    camera_t = 1:end_queue_row - start_queue_row +1 ;
    
    for track_row_num = start_queue_row:end_queue_row
        
        %% prepare image for warping, (were not frame for warping, directly write into final video sequence)
        % keep src_frame = fc frame
        if isempty(find(src_frame_next,1)) % first time read frame, read into src_frame_next
            if hasFrame(video_seq)
                src_frame_next = readFrame(video_seq);
                fc = fc + 1
                src_frame = src_frame_next;
                src_frame_next = readFrame(video_seq); % preload next frame
            end
        else
            if hasFrame(video_seq)
                src_frame_prev  =src_frame;
                src_frame = src_frame_next;
                src_frame_next = readFrame(video_seq);
                fc = fc + 1
            end
        end
        loop_flag = (fc ~= buf_queue(track_row_num).fc );
        while loop_flag == 1
            %imshow(src_frame);
            warp_img = double(src_frame);
            if fc <10
                write_prefix = 0;
            else
                write_prefix = [];
            end
            if WRITE_WARP == 1
                %imwrite(src_frame,[Folder '\f_' int2str(write_prefix) int2str(fc) '_nowarp.png']);
                writeVideo(Hyperlapse,src_frame);
            end
            if hasFrame(video_seq)
                src_frame_prev  =src_frame;
                src_frame = src_frame_next;
                src_frame_next = readFrame(video_seq);
                fc = fc + 1
            end
            loop_flag = (fc ~= buf_queue(track_row_num).fc );
        end
        
        %% prepare optimal warping homography matrix
        % record opt & orig camera path
        P_mat = P_queue(:,:,track_row_num);
        opt_x = [opt_x P_mat(1,3)];
        opt_y = [opt_y P_mat(2,3)];
        C_mat = buf_queue(track_row_num).F;
        camera_x = [camera_x C_mat(1,3)];
        camera_y = [camera_y C_mat(2,3)];
        
        dst = C_mat * src; % the dst used to warp back to 0
        C_inv = homography_4pts(dst(1:2,:), src(1:2,:));
        C_inv = C_inv./C_inv(3,3);
        B = P_mat * C_inv;
        H_inv = inv(B);
        H_inv = H_inv ./ H_inv(3,3);
        dst = B * src; % really target of warping
        
        minx = min(dst(1,:));   maxx = max(dst(1,:));
        miny = min(dst(2,:));   maxy = max(dst(2,:));
        
        T_prev2curr = F_queue(track_row_num).F; % prev frame is certainly in the same track
        if track_row_num+1 <= end_queue_row
            T_next2curr = inv(F_queue(track_row_num+1).F);
        else % case that next frame not in same track
            T_next2curr = eye(3);
        end
        
        %% real warping and write out to frame & seq
        if track_row_num+1 <= end_queue_row
            warp_img = stitching_3(src_frame_prev,src_frame,src_frame_next,T_prev2curr,T_next2curr,B);
        else % case that next frame not in same track
            warp_img = stitching_3(src_frame_prev,src_frame,0,T_prev2curr,T_next2curr,B);
        end

        if fc <10
            write_prefix = 0;
        else
            write_prefix = [];
        end
        if WRITE_WARP == 1
            %imwrite(uint8(warp_img),[Folder '\f_' int2str(write_prefix) int2str(fc) '_warp.png']);
            writeVideo(Hyperlapse,uint8(warp_img));
        end
    end
    if DRAW_PATH == 1
        camera_t = camera_t + start_end_of_track(curr_track,1) -1;
        figure
        plot3(camera_x,camera_y,camera_t,'ro-'),grid on,xlabel('x offset'),ylabel('y offset'),zlabel('t')
        hold on,plot3(opt_x,opt_y,camera_t,'b*-'),view([0 90 0])
    end
end

close(Hyperlapse);



