function [zhat,zvar]= one_bit_quant(phat,pvar,nuw,Y)

vecet=sqrt(2*((nuw^2)+pvar));
vecu=pvar./vecet;

eta1_r=sign(real(Y)).*real(phat);
eta1_r=2*eta1_r./vecet;

eta1_i=sign(imag(Y)).*imag(phat);
eta1_i=2*eta1_i./vecet;

eta1_r=min(1000,eta1_r);
eta1_i=min(1000,eta1_i);

phirat_r=normpdf(eta1_r)./normcdf(eta1_r);
phirat_i=normpdf(eta1_i)./normcdf(eta1_i);

zhat_r=(sign(real(Y)).*vecu).*phirat_r;
zhat_i=(sign(imag(Y)).*vecu).*phirat_i;
zhat=phat+zhat_r+(1i*zhat_i);

pvar_r=(eta1_r.*phirat_r)+(phirat_r.^2);
pvar_i=(eta1_i.*phirat_i)+(phirat_i.^2);
zvar=-(vecu.^2).*(pvar_r+pvar_i);
zvar=pvar+zvar;
end