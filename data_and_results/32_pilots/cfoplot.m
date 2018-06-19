clear all;
clc;
x=[-10:5:20];
y1=[3.9924    3.8879    1.0337   -2.5183   -8.0823  -12.9653  -34.2174];
y2=[-0.0704    0.3117   -2.6385  -26.6214  -55.8440  -56.5606  -56.6084];

plot(x,y1,'bd-','LineWidth',2,'MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','y')
hold on;
plot(x,y2,'rs-','LineWidth',2,'MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','b')

xlab=xlabel('$$\mathrm{SNR}$$(dB)','Interpreter','Latex')
set(xlab,'FontSize',18);
ylab=ylabel('CFO MSE (dB)');
set(ylab,'FontSize',18);
set(gca,'fontsize',16);
h_legend=legend('$$N_p=16, \Delta f_c=187.5\mathrm{KHz}$$','$$N_p=32, \Delta f_c=281.25\mathrm{KHz}$$');
set(h_legend,'Interpreter','Latex','FontSize',19);

grid on;