close all, clc, clear

% ------ get input from user
prompt = 'What is the original video file? ';
filename = input(prompt,'s')
video_seq=VideoReader(filename);
prompt = 'What is the output video file? ';
out_filename = input(prompt,'s')

% ----- modes
WRITE_FRAME = 1;
% ----- parameters
speedup_rate = 10; % speed-up rate
w = 2*speedup_rate;
fwidth = video_seq.Width;
fheight = video_seq.Height;
num_frame = floor(video_seq.FrameRate*video_seq.Duration); % count from 1 
d = sqrt(fwidth^2 + fheight^2); % in unit [pixel]
tou_c = 0.1*d;
gamma = 0.5*d;
alpha = 0.8;
tou_a = 200;
tou_s = 200;
g = 4;
lambda_s = 200;
lambda_a = 80;

global offset_x;
global offset_y;
global f_count_f;

offset_x = [];
offset_y = [];
f_count_f = 1;
% ----- cost matrix
% Cr, C0, V_i only store consecutive 2 frame cost
% Cm is the same as in paper
Cr = zeros(1,num_frame-1);
C0 = zeros(1,num_frame-1);
Cm = zeros(num_frame,num_frame);
Cm(Cm==0) = NaN;
V_i = zeros(1,num_frame-1);
Cs = zeros(num_frame,num_frame);
Cs(Cs==0) = NaN;


%% first time read begining w frames
no_use_term=[];
curr_t=1;
% keep w frames at a moment
for i=1:w+1
    if hasFrame(video_seq)
        eval(['frame_' int2str(i) '=readFrame(video_seq);']);
    end
end
for i=curr_t:curr_t+w-1
    %curr_t = i;
    eval(['[C_0,C_r,mu_i] = calc_CrC0(frame_' int2str(i) ',frame_' int2str(i+1) ',1,tou_c,gamma);']);
    Cr(1,i) = C_r; % eg: curr_t=1 means Cr(1,2), curr_t=2 means Cr(2,3)
    C0(1,i) = C_0;
    V_i(1,i) = mu_i;
end
C0_approx = cumsum(C0(1,curr_t:curr_t+w-1)); 
% approximation pair of C0(t,t+1),C0(t,t+2),C0(t,t+3).. C0(t,t+w-1)
mask = (C0_approx <= 0.05*d);
mask_c = (1-mask);
%{
    1.need to calc real Cr,C0 for pair(i,j) where apporx>0.05d
      only consider pair relating curr_t, namely 
      Cr(curr_t,curr_t+2),Cr(curr_t,curr_t+3),....,Cr(curr_t,curr_t+w-1)
    2.C0_no_approx save as same way as Cr, eg C0_no_approx(1,5) means C0
      among (1,6)
%}
no_approx_pos = find(mask_c);
C0_no_approx = zeros(1,w);
if ~isempty(no_approx_pos)
    no_approx_pos = no_approx_pos(1,1);
    for i=no_approx_pos:w
        eval(['[C_0,C_r,no_use_term] = calc_CrC0(frame_1,frame_' int2str(i+1) ',0,tou_c,gamma);']);
        if C_r<tou_c
            C0_no_approx(1,i) = C_0;
        else
            C0_no_approx(1,i) = gamma;
        end
    end
end
Cm(curr_t,curr_t+1:curr_t+w) = mask.*C0_approx + mask_c.*C0_no_approx; 
% no cost at Cm(t,t) 


%% following only read 1 new frame required, and only need calc 1 C_0,C_r,mu_i
%% this condition only invole all frame # that can calc w frame
while curr_t+w < num_frame
    curr_t = curr_t + 1
    end_t = curr_t + w;
    % shift frame & load 1 new in
    for i=1:w
        eval(['frame_' int2str(i) '=frame_' int2str(i+1) ';']);
    end
    if hasFrame(video_seq)
        eval(['frame_' int2str(w+1) '=readFrame(video_seq);']);
    else
        fprintf('no more frame could be read in! \n');
    end
    % calc new C_0 C_r mu_i for frame(end_t-1,end_t), newly computable
    % consecutive

    eval(['[C_0,C_r,mu_i] = calc_CrC0(frame_' int2str(w) ',frame_' int2str(w+1) ',1,tou_c,gamma);']);
    Cr(1,end_t-1) = C_r; % eg: curr_t=1 means Cr(1,2), curr_t=2 means Cr(2,3)
    C0(1,end_t-1) = C_0;
    V_i(1,end_t-1) = mu_i;
    
    C0_approx = cumsum(C0(1,curr_t:curr_t+w-1));
    % approximation pair of C0(t,t+1),C0(t,t+2),C0(t,t+3).. C0(t,t+w-1)
    mask = (C0_approx <= 0.05*d);
    mask_c = (1-mask);

    no_approx_pos = find(mask_c);
    C0_no_approx = zeros(1,w);
    if ~isempty(no_approx_pos)
        no_approx_pos = no_approx_pos(1,1);
        for i=no_approx_pos:w
            eval(['[C_0,C_r,no_use_term] = calc_CrC0(frame_1,frame_' int2str(i+1) ',0,tou_c,gamma);']);
            if C_r<tou_c
                C0_no_approx(1,i) = C_0;
            else
                C0_no_approx(1,i) = gamma;
            end
        end
    end
    Cm(curr_t,curr_t+1:curr_t+w) = mask.*C0_approx + mask_c.*C0_no_approx;

end


%% last few frames that cannot hold w frames for calculating, which means no more new
%% frame could be read
for minus_f = 1:w
    curr_t = curr_t + 1
    for i=1:w+1-minus_f
        eval(['frame_' int2str(i) '=frame_' int2str(i+1) ';']);
    end
    C0_approx = cumsum(C0(1,curr_t:end));
    % approximation pair of C0(t,t+1),C0(t,t+2),C0(t,t+3).. C0(t,t+w-1)
    mask = (C0_approx <= 0.05*d);
    mask_c = (1-mask);
    
    no_approx_pos = find(mask_c);
    C0_no_approx = zeros(1,w-minus_f);
    if ~isempty(no_approx_pos)
        no_approx_pos = no_approx_pos(1,1);
        for i=no_approx_pos:w-minus_f
            eval(['[C_0,C_r,no_use_term] = calc_CrC0(frame_1,frame_' int2str(i+1) ',0,tou_c,gamma);']);
            if C_r<tou_c
                C0_no_approx(1,i) = C_0;
            else
                C0_no_approx(1,i) = gamma;
            end
        end
    end
    Cm(curr_t,curr_t+1:end) = mask.*C0_approx + mask_c.*C0_no_approx;
end

%% interpolate V_i for NaN points & calc. Cs 
V_i(isnan(V_i)) = interp1(find(~isnan(V_i)), V_i(~isnan(V_i)), find(isnan(V_i)), 'cubic'); 
avg_camera_velocity = mean(V_i);
vi_in_Cs = alpha*speedup_rate*avg_camera_velocity./V_i + (1-alpha)*speedup_rate;
vi_in_Cs(1,:) = speedup_rate;
Cs = calc_Cs(num_frame,vi_in_Cs,tou_s,w);

%% DTW 
Cm = Cm*50;
[Dv,Dv_non_cum,Tv,path]=DTW(Cm,Cs,g,lambda_s,lambda_a,num_frame,tou_a,w);

%% write out selected frames
video_seq=VideoReader(filename); % reload src seq


if WRITE_FRAME == 1
    selectedFrame = VideoWriter(out_filename,'MPEG-4');
    open(selectedFrame);
    path_count = 1;
    f_count = 1;
    while f_count <= num_frame
        if hasFrame(video_seq)
            frame_buf =readFrame(video_seq);
        end
        if path(1,path_count)==f_count
            writeVideo(selectedFrame,frame_buf);
            path_count = path_count+1;
        end
        if path_count > numel(path)
            break;
        end
        f_count = f_count+1;
    end
    close(selectedFrame);
end

