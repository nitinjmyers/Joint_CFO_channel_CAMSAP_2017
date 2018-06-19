function Rate=ratecompute(H,Hest,SNR)
    [Uest,Sest,Vest]=svd(Hest);
    et=0.3634;
    npow=diag(Sest).*diag(Sest);
    npow=1./npow;
    scl=10^(-SNR/10);
    npow=npow*scl;
    Pt=1;%10^(SNR/10);
    palloc=wfill(npow.',Pt);
    Rs=diag(palloc);

    Nr=size(H,1);
    
    EfNsC=(scl*eye(Nr))+(et*diag(diag(H*Rs*H')));
    EfNsC=(1-et)*EfNsC;

    Pmat=Uest'*H*Vest;
    rk=size(Uest,2);
    Rate=0;
    for i=1:1:rk
     pow= (abs(Pmat(i,:)).^2)*((1-et)^2);
     pow=pow.*palloc;
     nsterm=abs(Uest(:,i)'*EfNsC*Uest(:,i))^2;
     cursnr=pow(i)/(nsterm+sum(pow)-pow(i));
     Rate=Rate+log2(1+cursnr);
    end
end

