clear all;
%clc;
nargin=0;
% To illustrate the performance of our algorithm for  a single run 
addpath('../main');
addpath('../EMGMAMP');
%% Test Case
Niter=1;
tic;
SNRvec=0:5:0;
Npil=64;
MSE=zeros(Niter,length(SNRvec));  % saves channel mean squared error
for ss=1:1:length(SNRvec)
    
CFOmat=zeros(Npil,Niter);   % saves b, i.e., the vector recovered after decomposing the lifted vector
for tt=1:1:Niter 
    
rng(21*tt); 
if nargin == 0
    
    %Fresh slate
    clc
    
    %Handle random seed
    if 1 %change to zero to try the same random draw repeatedly
        savedState = rng;
       % save random_state.mat savedState;
    else
        %load random_state.mat %#ok<UNRCH>
    end
    rng(savedState);
    
    %Control algorithms
    optIn.tryPbigamp = 0;
    optIn.tryEMPbigamp = 1;
    optIn.trySparseLift = 0;
    Nr=16;
    Nt=16;
    Ncl=2;
  

    [A,cfos,chs,ynf]=probsetup(Nr,Nt,Ncl,Npil); % to each measurement add sig*CN(0,1) later
    %Specify a problem size
    optIn.Nb = 1; %size of call vector
    optIn.Nc = Nr*Nt*Npil; %length of signal
    optIn.K = Ncl*Npil; %number of non-zeros in c
    optIn.M = Nr*Npil; %number of measurements
    
    flag=0;    
    %Specify SNR
    optIn.SNR = SNRvec(ss);
    
    %Control uniformVariance
    optIn.uniformVariance = 0;
    
else
    
    %Don't run stupid cases
    if optIn.M < (optIn.Nb + optIn.K)
        error('No reason to waste computer time')
    end
    
end


%% Problem setup

%Algorithm settings
tryPbigamp = optIn.tryPbigamp;
tryEMPbigamp = optIn.tryEMPbigamp;
trySparseLift = optIn.trySparseLift;

%Get params
Nb = optIn.Nb; %size of call vector
Nc = optIn.Nc; %length of signal
K = optIn.K; %number of non-zeros in c
M = optIn.M; %number of measurements
SNR = optIn.SNR;
%Determine nuw
nuw=sqrt(Nt)/sqrt(10^(SNR/10)); 
%Fix L
L = 1;


%% Build true signals
b = 1;
% The c below is the lifted vector "x" in the paper
c = vec(cfos*chs.');

%% Build the Z object

HO=ones(M,1);
%Build the A operator
Aop = Calibration_PLinTrans(A,HO);
zObject = Multiple_Snapshot_ParametricZ(Aop,L);


%% Continue setup
%Save the true Z
z = zObject.computeZ(b,c);

y = z + sqrt(nuw/2)*complex(randn(size(z)),randn(size(z)));
y=sign(real(y))+1i*sign(imag(y));
%% Establish the channel objects for P-BiG-AMP

%Prior on B
gB = CAwgnEstimIn(0, 1);
%Prior on C
gC = CAwgnEstimIn(0, 1);
gC = SparseScaEstim(gC,K/Nc);

%Output log likelihood
gOut = CAwgnEstimOut(y, nuw);


%% Options and Problem objects

%Setup the problem
problem = PBiGAMPProblem();
problem.M = M;
problem.Nb = Nb;
problem.Nc = Nc;
problem.zObject = zObject;

%Setup the options
opt = PBiGAMPOpt();

%Uniform variance
opt.uniformVariance = optIn.uniformVariance;

%Iteration settings
opt.verbose =1;% (nargin == 0);
opt.nit = 30;
opt.tol = 1e-4;

%Specify maximum number of trials
maxTrials = 1;     % was 3 for CAMSAP paper
threshNMSE = max(-15-SNR,-60);

%Use the initial values
opt.bhat0 = 1;%sqrt(1/2)*(randn(Nb,1) + 1j*randn(Nb,1));
opt.chat0 = sqrt(1/2)*(randn(Nc,1) + 1j*randn(Nc,1));

%Variances
[~,holder,~] = gB.estimInit();
opt.bvar0 = 10*holder;
[~,holder,~] = gC.estimInit();
opt.cvar0 = 10*holder;

%Specify error functions
opt.error_functionB = @(qval) 0;
opt.error_functionC = @(qval) 0;
opt.error_function = @(qval) 20*log10(norm(qval - y,'fro') / norm(y,'fro'));

%Overall error function
error_function_full = @(bhat,chat) 20*log10(norm(bhat*chat.' - b*c.','fro')/norm(b*c.','fro'));


%% Init results

%Initialize results as empty
results = [];



%% P-BiG-AMP

if tryPbigamp
    
    %Count failures
    failCounter = 0;
    stop = false;
    
    tstart = tic;
    while ~stop
        
        %Increment counter
        failCounter = failCounter + 1; 
        [estFin, ~, estHist] = ...
            PBiGAMP(gB, gC, gOut, problem, opt);
        
        
        if (estHist.errZ(end) < threshNMSE) || failCounter > maxTrials
            stop = true;
        else
            opt.bhat0 = 1;%sqrt(1/2)*(randn(Nb,1) + 1j*randn(Nb,1));
            opt.chat0 = sqrt(1/2)*(randn(Nc,1) + 1j*randn(Nc,1));
        end
        disp(['Attempts completed: ' num2str(failCounter) ' Final Z error= ' ...
            num2str(estHist.errZ(end))])
       
    end
    tGAMP = toc(tstart);
    
    
    %Save results
    loc = length(results) + 1;
    results{loc}.name = 'P-BiG-AMP'; %#ok<*AGROW>
    results{loc}.err = estHist.errZ(end);
    results{loc}.errFull = error_function_full(estFin.bhat,estFin.chat);
    results{loc}.time = tGAMP;
    results{loc}.errHist = estHist.errZ;
    results{loc}.timeHist = estHist.timing;
    results{loc}.fails = failCounter;
end

%% Try EM-P-BiG-AMP

if tryEMPbigamp
    
    %Specify options
    EMopt.B_type = 'CG';
    EMopt.C_type = 'CBG';
    EMopt.C_var_init = 'fixed';
    EMopt.C_learn_lambda = 1;
    EMopt.C_learn_var = 1;
    EMopt.learn_noisevar=0;
    EMopt.C_lambda = 0.01; %start smaller
 %   EMopt.B_lambda = 0.01; %start smaller
    
    EMopt.maxEMiter = 1;  % was 4 in CAMSAP paper
    opt.nit = 30;
    EMopt.noise_var=nuw;
    %Coding
%    disp('Starting EM-P-BiG-AMP')
    olderr=100;
    %Count failures
    failCounter = 0;
    stop = false;
    
    tstart = tic;
    while ~stop
        
        %Increment counter
        failCounter = failCounter + 1;
        
        %Run EM-P-BiG-AMP
        [estFinEM, ~, ~,estHistEM] = ...
            EMPBiGAMP(y,problem,opt,EMopt);
       
        if(estHistEM.errZ(end)<olderr)   %can do this -checks norm(y-Ax)
                flag=1;
                olderr=estHistEM.errZ(end);
                cold=estFinEM.chat;
                bold=estFinEM.bhat;
        end
        %Check
        if (estHistEM.errZ(end) < threshNMSE) || failCounter > maxTrials
            stop = true;
            
        else      
            opt.bhat0 = sqrt(1/2)*(randn(Nb,1) + 1j*randn(Nb,1));
            opt.chat0 = sqrt(1/2)*(randn(Nc,1) + 1j*randn(Nc,1));
        end
        gest=estFinEM.chat;
        disp(['Attempts completed: ' num2str(failCounter) ' Final Z error= ' ...
            num2str(estHistEM.errZ(end))])
        flag
     if(flag==1)
           stop=true;
     end
    end
    tEMGAMP = toc(tstart);
    
    loc = length(results) + 1;
    results{loc}.name = 'EM-P-BiG-AMP'; %#ok<*AGROW>
    results{loc}.err = estHistEM.errZ(end);
    results{loc}.errFull = error_function_full(estFinEM.bhat,estFinEM.chat);
    results{loc}.time = tEMGAMP;
    results{loc}.errHist = estHistEM.errZ;
    results{loc}.timeHist = estHistEM.timing;
    results{loc}.fails = failCounter;
end

if(flag==1)
    gest=cold;
else
    gest=zeros(Npil*Nr*Nt,1);
end
Cest=reshape(gest,[Npil,Nr*Nt]);   % this is an estimate of the lifted vector
[U,S,V]=svd(Cest);
CFOmat(:,tt)=U(:,1);
CHSP=conj(V(:,1));

sc=CHSP'*chs/(norm(CHSP,'fro')^2);
MSE(tt,ss)=(norm(chs-(sc*CHSP),'fro')/norm(chs,'fro'))^2;
CHmat=dftmtx(Nr)*reshape(chs,[Nr,Nt])*dftmtx(Nt)'/sqrt(Nr*Nt);
CHest=dftmtx(Nr)*reshape(sc*CHSP,[Nr,Nt])*dftmtx(Nt)'/sqrt(Nr*Nt);
end
%dlmwrite(strcat('final_64pilotsCFO_vec_SNR_',num2str(SNRvec(ss)),'.txt'),CFOmat);
end
toc
figure(1)
subplot(1,2,1)
surf(abs(reshape(chs,[Nr,Nt])))
view([0,90])
title('True beamspace channel')
subplot(1,2,2)
surf(abs(reshape(sc*CHSP,[Nr,Nt])))
view([0,90])
title('Estimated beamspace channel')

figure(2)
plot(abs(CFOmat(:,1)))
title('Estimated DFT of phase error vector due to CFO')
sprintf('Channel NMSE(dB)=%f',10*log10(MSE))

%dlmwrite('final_64pilots_MSE.txt',MSE)