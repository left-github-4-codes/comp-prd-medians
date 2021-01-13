hold off;
x=H(1,1:5);
%y=-1.7317456-0.8184845.*x; %six points deepest line
y=-2.087114-1.009444.*x; %five points deepest line
plot(x,y, 'LineWidth', 1,'DisplayName', 'Line L1', 'color', 'r');
hold on;
%x=[-7 12];
y=-1.87-0.977.*x;
plot(x,y, 'LineWidth', 1, 'DisplayName', 'Line L2', 'color', 'b');
%x=[-7 12];
%y=0.07-0.08.*x;
y=-1.58-0.77.*x; %deepest RD line 3/5
%disp([x,y]);
plot(x,y, 'LineWidth', 1,'DisplayName', 'Line L3','color', 'g');
%legend({'Line L1','Line L2', 'Line L3'},'Location','southwest');
%scatter([-7 12], [-5 5]);
%H1=H(1:5,:);
scatter(H(1,1:5),H(2,1:5),'filled');
lgd=legend({'Line L1','Line L2', 'Line L3'},'Location','southwest');
title(lgd, 'Three regression lines')

