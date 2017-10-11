clc

prompt = 'What is the video file name to be cropped? ';
%seq_to_crop = input(prompt,'s')
seq_to_crop = 'hyperlapse_bike_Cm50_sw75.mp4';
video_seq=VideoReader(seq_to_crop);
prompt = 'What is the final video file name? ';
final_name = input(prompt,'s')
final_seq = VideoWriter(final_name,'MPEG-4');

prompt = 'up crop ratio [%]? ';
crop_up_ratio = input(prompt);
prompt = 'bottom crop ratio [%]? ';
crop_bottom_ratio = input(prompt);
prompt = 'left crop ratio [%]? ';
crop_left_ratio = input(prompt);
prompt = 'right crop ratio [%]? ';
crop_right_ratio = input(prompt);

Folder = 'Crop_test';

%
crop_up_pel = floor(video_seq.Height * crop_up_ratio/100);
crop_bottom_pel = floor(video_seq.Height * crop_bottom_ratio/100);
crop_left_pel = floor(video_seq.Width * crop_left_ratio/100);
crop_right_pel = floor(video_seq.Width * crop_right_ratio/100);

final_seq.FrameRate = 30;
open(final_seq);

fc = 0;
while fc < video_seq.FrameRate * video_seq.Duration
    if hasFrame(video_seq)
        orig_frame = readFrame(video_seq);
        fc = fc+1
    else
        break
    end
    cropped_frame = orig_frame(crop_up_pel:end-crop_bottom_pel, crop_left_pel:end-crop_right_pel, :);
    writeVideo(final_seq,cropped_frame);
    imwrite(cropped_frame,[Folder '\f_' int2str(fc) '.png']);
end

close(final_seq);