display_num = 100;
[y,x] = meshgrid(g:display_num,g:display_num);
surf(x,y,Dv_non_cum(1:display_num-g+1,1:display_num-g+1)),colormap(jet)
xlabel('current'),ylabel('next'),zlabel('cost')
view([90 90])
hold on

xl = path(1,1:end-1);
xl(xl>display_num) = [];
yl = path(1,2:end);
yl(yl>display_num) = [];
if numel(xl) > numel(yl)
    xl = xl(1,1:end-1);
end
tmp = ones(1,numel(xl));
plot3(xl,yl,tmp,'ro-','MarkerSize',10,'LineWidth',2);

set(gca,'FontSize',25)