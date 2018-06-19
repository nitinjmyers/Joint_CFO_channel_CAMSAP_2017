clear all;
clc;

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
semilogy(mean((W-we).*(W-we)),'ro-')
%plot(10*log10(mean((W-we).*(W-we))),'ro-')

