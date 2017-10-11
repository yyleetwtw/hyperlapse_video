mask = (Cr==100000);
no_T_pos = find(mask);

padded_offset_x = zeros(1,numel(Cr));
padded_offset_y = zeros(1,numel(Cr));

for i = 1:numel(no_T_pos)
    insert_pos = no_T_pos(1,i);
    if i==1
        padded_offset_x(1,1:insert_pos-1) = offset_x(1,1:insert_pos-1);
        padded_offset_x(1,insert_pos) = 100000;
        padded_offset_y(1,1:insert_pos-1) = offset_y(1,1:insert_pos-1);
        padded_offset_y(1,insert_pos) = 100000;
        end_pos_in_orig_offset = insert_pos; % yet used data
    else
        last_padded_pos = no_T_pos(1,i-1) +1; % yet used data
        fill_in_num = insert_pos - last_padded_pos ;
        last_end_pos_in_orig_offset = end_pos_in_orig_offset;
        end_pos_in_orig_offset = end_pos_in_orig_offset + fill_in_num;
        if last_end_pos_in_orig_offset ~= end_pos_in_orig_offset % means fill in required
            padded_offset_x(1,last_padded_pos:last_padded_pos+fill_in_num-1) = ...
                offset_x(1,last_end_pos_in_orig_offset:end_pos_in_orig_offset-1);
            padded_offset_y(1,last_padded_pos:last_padded_pos+fill_in_num-1) = ...
                offset_y(1,last_end_pos_in_orig_offset:end_pos_in_orig_offset-1);
        end
        padded_offset_x(1,insert_pos) = 100000;
        padded_offset_y(1,insert_pos) = 100000;
    end
end
padded_offset_x(1,no_T_pos(1,end)+1:end) = offset_x(1,end_pos_in_orig_offset:end);
padded_offset_y(1,no_T_pos(1,end)+1:end) = offset_y(1,end_pos_in_orig_offset:end);

offset_sum = sqrt(padded_offset_x.^2+padded_offset_y.^2);
figure,plot(1:10:numel(Cr),offset_sum(1,1:10:numel(Cr)),'k-','LineWidth',4),hold on
plot(path,offset_sum(path+1),'r-','LineWidth',4)
grid on, xlabel('frame count'),ylabel('offset'),legend('naive skip','selected path')
set(gca,'FontSize',30)


%{
figure, plot(path,padded_offset_x(path+1),'ro-'),hold on
plot(1:10:numel(Cr),padded_offset_x(1,1:10:numel(Cr)),'k*-')
grid on, xlabel('frame count'),ylabel('x offset'),legend('selected','naive')
figure, plot(path,padded_offset_y(path+1),'ro-'),hold on
plot(1:10:numel(Cr),padded_offset_y(1,1:10:numel(Cr)),'k*-')
grid on, xlabel('frame count'),ylabel('y offset'),legend('selected','naive')
%}