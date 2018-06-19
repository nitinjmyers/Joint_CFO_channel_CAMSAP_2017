clear all;
clc;
Nr=32;
Nt=32;
Ncl=2;

Np=64;
sf=1;
qb=3;  %quantization bits
pn_sig=0;
sig=0.5;
Gr=sf*Nr;
Gt=sf*Nt;
UNt=ntnmtx(Gt,Nt);
UNr=ntnmtx(Gr,Nr);
H=Hmat(Nr,Nt,Ncl);
% Training sequence
T=zeros(Nt,Np);

T= randi(2^qb,[Nt,Np])*2*pi/(2^qb);
T=exp(1i*T)/sqrt(Nt);
f_s=2000;   %KHz
f_e=31.25*8.5;
we=2*pi*f_e/f_s;
cfov=exp(1i*we*[0:1:Np-1]).';
%fe_max=280;  %KHz less than fs/2 - 40 meas considered
%gamma=3;
%wmax=2*pi*fe_max*gamma/f_s  ; % from CFO limits
%Nmax=floor(Np*wmax/(2*pi));
B=dftmtx(Np)'/sqrt(Np);  % basis for CFO
%B=[B(:,1:Nmax),B(:,end:-1:Np-Nmax+1)];   % from knowledge of worst case CFO 0-->f_e
kc=size(B,2);
pnv=zeros(Np,1);
pnv(1)=pn_sig*randn;
for i=2:1:Np
    pnv(i)=pnv(i-1)+ pn_sig*randn;
end 
pnv=exp(1i*pnv);
noise=(randn(Nr,Np)+1i*randn(Nr,Np))/sqrt(2);
Y=H*T*diag(cfov.*pnv)+(sig*noise);
Y=sign(real(Y))+1i*sign(imag(Y));
A=kron(UNr,(UNt'*T).');
y=vec(Y.');
Bf=kron(ones(Nr,1),B);
Am=zeros(Np*Nr,kc*Gt*Gr);
for i=1:1:Np*Nr
    Am(i,:)=kron(A(i,:),Bf(i,:));
end
ho=B'*(cfov.*pnv);
xo=vec((UNr'*H*UNt).');
fo=ho*xo.';
norm(y-Am*vec(fo))
tic;
xetp=solve_OMP(y,Am,sig*sqrt(Nr*Np),300);
toc
fest=reshape(xetp,[kc,Gt*Gr]);
[a,b,c]=svd(fest);
Xest=reshape(conj(c(:,1)),[Gt,Gr]).';
surf(abs(Xest))
view([0 90]);
figure()
surf(abs(UNr'*H*UNt))
view([0 90]);

