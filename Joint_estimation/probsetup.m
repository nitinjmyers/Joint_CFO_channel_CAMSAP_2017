function [Am,ho,xo,y]=probsetup(Nr,Nt,Ncl,Npil)
pn_sig=0;
UNt=dftmtx(Nt)/sqrt(Nt);
UNr=dftmtx(Nr)/sqrt(Nr);
H=Hmat(Nr,Nt,Ncl);
% Training sequence
T=zeros(Nt,Npil);
qb=2;
T= randi(2^qb,[Nt,Npil])*2*pi/(2^qb);
T=exp(1i*T)*exp(pi/2^qb);
%T=(randn(Nt,Npil)+1i*randn(Nt,Npil))/sqrt(2);
f_s=2000;   %KHz
f_e=2000*3.5/Npil;
we=2*pi*f_e/f_s;
cfov=exp(1i*we*[0:1:Npil-1]).';
%fe_max=280;  %KHz less than fs/2 - 40 meas considered
%gamma=3;
%wmax=2*pi*fe_max*gamma/f_s  ; % from CFO limits
%Nmax=floor(Npil*wmax/(2*pi));
B=dftmtx(Npil)'/sqrt(Npil);  % basis for CFO
%B=[B(:,1:Nmax),B(:,end:-1:Npil-Nmax+1)];   % from knowledge of worst case CFO 0-->f_e
kc=size(B,2);
pnv=zeros(Npil,1);
pnv(1)=pn_sig*randn;
for i=2:1:Npil
    pnv(i)=pnv(i-1)+ pn_sig*randn;
end 
pnv=exp(1i*pnv);
%noise=(randn(Nr,Npil)+1i*randn(Nr,Npil))/sqrt(2);
Y=H*T*diag(cfov.*pnv);
%Y=sign(real(Y))+1i*sign(imag(Y));
A=kron(UNr,(UNt'*T).');
y=vec(Y.');
Bf=kron(ones(Nr,1),B);
Am=zeros(Npil*Nr,kc*Nt*Nr);
for i=1:1:Npil*Nr
    Am(i,:)=kron(A(i,:),Bf(i,:));
end
ho=B'*(cfov.*pnv);
xo=vec((UNr'*H*UNt).');
fo=ho*xo.';


end