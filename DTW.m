
function [Dv,Dv_non_cum,Tv,path]=DTW(Cm,Cs,g,lambda_s,lambda_a,num_frame,tou_a,w)

Dv = zeros(num_frame,num_frame);
Dv(Dv==0) = NaN;
Dv_non_cum = zeros(num_frame,num_frame);
Dv_non_cum(Dv_non_cum==0) = NaN;
Tv = zeros(num_frame,num_frame);
Tv(Tv==0) = NaN;
path = [];

% construct possible Ca
tmp=1:w;
j_minus_i=repmat(tmp,[w 1]);
tmp=w:-1:1;
i_minus_h=repmat(tmp',[1 w]);
tmp = (j_minus_i - i_minus_h).^2;
Ca = min(tmp,tou_a);

% initialization
Dv(1:g,:) = Cm(1:g,:) + lambda_s*Cs(1:g,:);

% first pass , row-based design
for f_num = g:1:num_frame-1
    horizontal_size = min(num_frame-f_num,w);
    c = Cm(f_num,f_num+1:f_num+horizontal_size) + lambda_s*Cs(f_num,f_num+1:f_num+horizontal_size);
    get_rows = min(f_num-1,w);
    tmp = Dv(f_num-1:-1:f_num-get_rows,f_num);
    tmp = repmat(tmp,[1 horizontal_size]);
    tmp = tmp + lambda_a*Ca(end:-1:end-get_rows+1,1:horizontal_size);
    [min_D,min_k] = min(tmp,[],1); % take min along col direction
    Dv(f_num,f_num+1:f_num+horizontal_size) = min_D+c;
    Tv(f_num,f_num+1:f_num+horizontal_size) = f_num-min_k;
    for non_cum_pre = 1:horizontal_size
        Dv_non_cum(f_num,f_num+non_cum_pre) = Dv(f_num,f_num+non_cum_pre) - ...
                                              Dv(Tv(f_num,f_num+non_cum_pre),f_num);
    end
end

% second pass, trace back
for i = num_frame-g:num_frame-1
    for j = 1:w
        if i == num_frame-g && j == 1
            min_D = Dv(i,i+j);
            s = i;
            d= i+j;
            continue
        end
        if i+j <=50
            tmp = Dv(i,i+j);
            if tmp<min_D
                min_D = tmp;
                s = i;
                d = i+j;
            end
        end
    end
end
path = [path d];
while s > g
    path = [path s];
    b = Tv(s,d);
    d = s;
    s = b;
end
path = [path s];
path = sort(path);

% normalize cost of Dv_non_cum at each row
for f_num = g:1:num_frame-1
    mask = (~isnan(Dv_non_cum(f_num,:)));
    curr_row = Dv_non_cum(f_num,:);
    row_max = max(curr_row(mask));
    curr_row = curr_row / row_max;
    Dv_non_cum(f_num,:) = curr_row;
end
    