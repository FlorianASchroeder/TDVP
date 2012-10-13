function [modelpara]=SBM_genpara(modelpara)

%Calculate the wilson chain parameters for spin-boson model
%See Appendix A of PRB 71, 045122 for details.
z=modelpara.z;

%J=@(w) 2*pi*modelpara.alpha*w.^modelpara.s;
fac=2*pi*modelpara.alpha/(1+modelpara.s);
tempfac=1/((1+modelpara.s)*log(modelpara.Lambda));
tempLam=modelpara.Lambda^(1+modelpara.s);
tempexp=1/(1+modelpara.s);

bigL=floor(-1*log(realmin)/log(modelpara.Lambda))/2 %Use a large enough start H to make sure the accuracy after transformed to Wilson chain

xi=zeros(bigL,1);
gamma=zeros(bigL,1);

%w0=1;
%wdata=zeros(bigL+1,1);
%wdata(1)=1;
for j=2:bigL+1
    %w1=modelpara.Lambda.^(2-j-z);
    %wdata(j)=w1;
    %gamma(j-1)=sqrt(quad(J,w1,w0));
    if j==2
	    xi(j-1)=(tempfac*(1-tempLam^(-z))-z+1)^tempexp;
	    gamma(j-1)=sqrt(fac*(1-tempLam^(-z)));
    else
    	xi(j-1)=(tempfac*tempLam^(2-j-z)*(tempLam-1))^tempexp;
	gamma(j-1)=sqrt(fac*(tempLam-1)*tempLam^(2-j-z));
    end
    %w0=w1;
end

%fprintf('wdata(end-5:end) = %.10g\n', wdata(end-5:end));
%bar(log(wdata),ones(modelpara.L+1,1));
% gamma(end-5:end)
% xi(end-5:end)

indiag=zeros(bigL,1);
inrow=indiag;

for n=2:bigL
    inrow(n)=gamma(n-1)/(2*sqrt(pi));
    indiag(n)=xi(n-1);
end
%inrow
%indiag
[epsilon,t]=star2tridiag(indiag,inrow);
[epsilon,t]=extroplate(epsilon,t,modelpara.L);
modelpara.epsilon=epsilon(1:modelpara.L-1);
modelpara.t=t(1:modelpara.L-1);
%modelpara.t=abs(modelpara.t); %I noticed that the r form hess() function sometimes changed sign.
