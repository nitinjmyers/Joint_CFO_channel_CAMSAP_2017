clear all;
clc;
plotflag=1; % use 1 ,2 and 3 for different plots
cd ../
load('500_iter_64_meas_1bit_CHP_MSqE_ignoreCFO.mat')
MSEign=MSE;
Rateign=Rate;
load('500_iter_64_meas_1bit_CHP_MSqE_noCFO.mat')
MSE64=MSE;%
Rate64=[Rate(:,1:6),Rate(:,6)];;
cd plots_includes_64_pilots
SNR=[-10:5:20];
SNRL=length(SNR);
sets=[40,50,70,150,200]; %number of iterations

Npil=64;
Nit=510;
CFO=[];
FullRate=[];
FullMSE=[];
for i=1:1:length(sets)
    load(strcat(num2str(sets(i)),'_iter_64_meas_1bit_CHP_MSqE_CFO.mat'));
    CFO=[CFO;CFOmat];
    FullRate=[FullRate;Rate];
    FullMSE=[FullMSE;MSE];
end

for jj=1:1:SNRL
    for ii=1:1:Nit
    N=2*Npil;
    sp=vec(CFO(ii,jj,:));
    sig=ifft(sp);
    sp=fft(sig,N);
    [val,kp]=max(abs(sp));
    lt=kp-1;
    rt=kp+1;
    if(lt==0)
        lt=N;
    end
    if(rt==N+1)
        rt=1;
    end
    corr=(tan(pi/N))/(pi/N);
    Dr=(2*sp(kp))-(sp(lt)+sp(rt));
    corr=corr*real((sp(lt)-sp(rt))/Dr);
    west=2*pi*(kp-1)/N;
    binwid=2*pi/N;
  
    west=west+(binwid*corr);
    WEM(jj,ii)=west;
    end
end
W=WEM.';
we=2*pi*3.5/Npil; % 1.5 or 3.5 depending on Npil- set according to paper - maximally off grid

W(find(W>pi))=W(find(W>pi))-(2*pi);
switch plotflag
    case 0
    h1=semilogy(SNR,mean((W-we).*(W-we)),'bd-','LineWidth',2,'MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','y')
    ylab=ylabel('MSE of CFO estimate');
    case 1
    h1=plot(SNR,10*log10(mean((W-we).*(W-we))),'bd-','LineWidth',2,'MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','y')
    ylab=ylabel('MSE of CFO estimate (dB)');
    case 2
    h1=plot(SNR,mean(FullRate),'bd-','LineWidth',2,'MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','y')
    ylab=ylabel('Achievable Rate (bps/Hz)');
    case 3
    h1=plot(SNR,10*log10(mean(FullMSE)),'bd-','LineWidth',2,'MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','y')
    ylab=ylabel('MSE of channel estimate (dB)');
end
hold on;
xlab=xlabel('$$\mathrm{SNR}$$(dB)','Interpreter','Latex')
set(xlab,'FontSize',18);
addpath('../32_pilots');

SNR=[-10:5:20];
SNRL=length(SNR);
sets=[200,500]; %number of iterations
Npil=32;
Nit=700;
CFO=[];
FullRate=[];
FullMSE=[];
for i=1:1:length(sets)
    load(strcat(num2str(sets(i)),'_iter_32_meas_1bit_CHP_MSqE_CFO.mat'));
    CFO=[CFO;CFOmat];
    FullRate=[FullRate;Rate];
    FullMSE=[FullMSE;MSE];
end

for jj=1:1:SNRL
    for ii=1:1:Nit
    N=1*Npil;
    sp=vec(CFO(ii,jj,:));
    sig=ifft(sp);
    sp=fft(sig,N);
    [val,kp]=max(abs(sp));
    lt=kp-1;
    rt=kp+1;
    if(lt==0)
        lt=N;
    end
    if(rt==N+1)
        rt=1;
    end
    corr=(tan(pi/N))/(pi/N);
    Dr=(2*sp(kp))-(sp(lt)+sp(rt));
    corr=corr*real((sp(lt)-sp(rt))/Dr);
    west=2*pi*(kp-1)/N;
    binwid=2*pi/N;
  
    west=west+(binwid*corr);
    WEM(jj,ii)=west;
    end
end
W=WEM.';
we=2*pi*1.5/Npil; % 1.5 or 3.5
W(find(W>pi))=W(find(W>pi))-(2*pi);
mean(Rate)
switch plotflag
    case 0
    h2=semilogy(SNR,mean((W-we).*(W-we)),'rs-','LineWidth',2,'MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','b')
    case 1
    h2=plot(SNR,10*log10(mean((W-we).*(W-we))),'rs-','LineWidth',2,'MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','b')
    case 2
    h2=plot(SNR,mean(FullRate),'rs-','LineWidth',2,'MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','b')
    case 3
    h2=plot(SNR,10*log10(mean(FullMSE)),'rs-','LineWidth',2,'MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','b')
end
set(ylab,'FontSize',18);
set(gca,'fontsize',16);

switch plotflag
    case 2
    h4=plot(SNR,mean(Rateign),'kx-','LineWidth',2,'MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','b')
    h3=plot(SNR,mean(Rate64),'mo--','LineWidth',2,'MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','g')      
    case 3
    h4=plot(SNR,10*log10(mean(MSEign)),'kx-','LineWidth',2,'MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','b')
    h3=plot(SNR,10*log10(mean(MSE64)),'mo--','LineWidth',2,'MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','g')
end

switch plotflag
    case 0
    h_legend=legend([h2 h1],'$$N_p=32, \Delta f_c=93.75 \,\mathrm{KHz}$$','$$N_p=64, \Delta f_c=109.375\,\mathrm{KHz}$$');
    case 1
    h_legend=legend([h2 h1],'$$N_p=32, \Delta f_c=93.75 \,\mathrm{KHz}$$','$$N_p=64, \Delta f_c=109.375\,\mathrm{KHz}$$');
    case 2
    h_legend=legend([h3 h1 h2 h4],'Oracle CFO, $$N_p=64$$','Proposed method, $$N_p=64$$', 'Proposed method, $$N_p=32$$','GAMP ignoring CFO, $$N_p=64$$' );
    case 3
    h_legend=legend([h4 h2 h1 h3],'GAMP ignoring CFO, $$N_p=64$$', 'Proposed method, $$N_p=32$$', 'Proposed method, $$N_p=64$$','Oracle CFO, $$N_p=64$$');
end
%ylim([-10,0.5])
set(h_legend,'Interpreter','Latex','FontSize',18);
