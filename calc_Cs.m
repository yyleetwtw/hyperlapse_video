function Cs = calc_Cs(num_frame,vi_in_Cs,tou_s,w)

Cs = zeros(num_frame,num_frame);
Cs(Cs==0) = NaN;
tmp=1:num_frame;
i=repmat(tmp',[1 num_frame]);
j=repmat(tmp,[num_frame 1]);
tmp=j-i;
vi_padded = repmat(vi_in_Cs',[1 num_frame]);
pader = NaN;
pader = repmat(pader,[1 num_frame]); % add last raw that not exist
vi_padded = [ vi_padded ; pader];
tmp = (tmp-vi_padded).^2;

% construc mask for Cs only (i,i+1:i+w-1)
mask = ones(num_frame,num_frame);
mask = triu(mask,1) - triu(mask,w+1);
a = min(tmp,tou_s);
Cs(logical(mask))=a(logical(mask));