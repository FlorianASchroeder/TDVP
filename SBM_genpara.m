function [modelpara]=SBM_genpara(modelpara)
%Calculate the wilson chain parameters for spin-boson model using logarithmic discretization.
%See Appendix A of Bulla 2005 PRB 71, 045122 for details.
%
% Modified:
%   FS 19/02/2014:  - added automatic L for L = 0;
%	FS 25/05/2014:	- added Spectral function of defined in Renger 2002,
%                     following Müh 2012
%	FS 27/05/2014:	- Changed to J(w) as in Renger 2002. This defines for
%                     B777, a monomer of B850, while Müh 2012 needs a
%                     prefactor to use it for pigments in PSII

z=modelpara.z;
bigL=ceil(floor(-1*log(realmin)/log(modelpara.Lambda))/2) %Use a large enough start H to make sure the accuracy after transformed to Wilson chain

function y = J_Renger(w,i)
    %% Section for Spectral function for LH2 of Rhodoblastus acidophilus from Renger 2002, Müh 2012
    % need for numerical integration.
    % multiplied by \pi * \omega^2: \pi to equal out division of \sqrt\pi later in the \gamma. Problem due to different
    %               definition of J(w) in Leggett et al.
    %               \omega^2 because Leggett defined J(w) as spectrum of the dimension-ful coupling constant! Renger
    %               took one \hbar\omega out.
    % This J(w) gives a Huang-Rhys factor of 1.3
    % i to get function J*w^(-i)

%    w=0:cmToeV(1):cmToeV(3000); i=0;                                                                % testing
    y = pi/(factorial(7)*2).*((0.8/cmToeV(0.56).^4 .* w.^(5-i).*exp(-(w/cmToeV(0.56)).^(1/2)))+(0.5/(cmToeV(1.94).^4) .* w.^(5-i).*exp(-(w/cmToeV(1.94)).^(1/2))));
%    plot(y);xlabel('$\omega / cm^{-1}$');ylabel('$J(\omega)/eV$'); formatPlot(1);                             % testing

    % original definition: J = @(w) 1/(factorial(7)*2).*((0.8/cmToeV(0.56).^4 .* w.^3.*exp(-(w/cmToeV(0.56)).^(1/2)))+(0.5/cmToeV(1.94).^4 .* w.^3.*exp(-(w/cmToeV(1.94)).^(1/2))));
    % from Cheng:
    % J=@(w) 2*pi*modelpara.alpha*w.^modelpara.s;
end

if strcmp(modelpara.model,'MLSpinBoson') && modelpara.MLSB_mode == 2
    %% use J(w) from: F. Müh and T. Renger, Biochim. Biophys. Acta 1817, 1446 (2012). DOI: 10.1016/j.bbabio.2012.02.016
    % this code section can be used with any definition of J(w)
    % Units: eV

%    tic; modelpara.Lambda = 2; bigL = 511; z=1;        % testing
%   w_c = 1000; accurate up to 10^-3 eV; w_c = 2000: 10^-6 eV; w_c = 3000: 10^-8 eV
    w_cutoff = cmToeV(3000);                            % upper cutoff frequency

    xi=zeros(bigL,1);
    gamma=zeros(bigL,1);
    w_limits = modelpara.Lambda.^(1-z-(0:bigL)).*w_cutoff;      % Define limits of all the discretization intervals
    w_limits(1) = w_cutoff;
    for j=1:bigL
        intJ = integral(@(w) J_Renger(w,0),w_limits(j+1),w_limits(j));
        intJoverW = integral(@(w) J_Renger(w,1),w_limits(j+1),w_limits(j));
        xi(j)= intJ/intJoverW;                          % bath energy levels
        gamma(j)= 2*sqrt(intJ);                         % bath couplings. Factor 2 because I don't want term \sigma_z/2. 2 is divided later again.
    end
    xi(isnan(xi))=0;                                    % if numbers < e-161 they become NaN. This is fix for it.
    gamma(isnan(gamma))=0;

%    toc(tic);                                          % testing
else
    %% use J(w) from A. Leggett et al., Rev. Mod. Phys. 59, 1 (1987). DOI: 10.1103/RevModPhys.59.1
    %           and R. Bulla, N.-H. Tong, and M. Vojta, Phys. Rev. Lett. 91, 170601 (2003). DOI: 10.1103/PhysRevLett.91.170601
    %           and R. Žitko and T. Pruschke, Phys. Rev. B 79, 085106 (2009). DOI: 10.1103/PhysRevB.79.085106
    % TODO: bring onto form as above, to make readable and shorter!
%    tic;modelpara.alpha = 0.5; modelpara.s = 1; modelpara.Lambda = 2; bigL = 511; z=1;       % testing
    fac=2*pi*modelpara.alpha/(1+modelpara.s);
    tempfac=1/((1+modelpara.s)*log(modelpara.Lambda));
    tempLam=modelpara.Lambda^(1+modelpara.s);
    tempexp=1/(1+modelpara.s);

    xi=zeros(bigL,1);
    gamma=zeros(bigL,1);

    % calculate discretized energy levels and couplings
    for j=2:bigL+1
        if j==2
            xi(j-1)=(tempfac*(1-tempLam^(-z))-z+1)^tempexp;
            gamma(j-1)=sqrt(fac*(1-tempLam^(-z)));                      % modify this after the vector calculation.
        else
            xi(j-1)=(tempfac*tempLam^(2-j-z)*(tempLam-1))^tempexp;
            gamma(j-1)=sqrt(fac*(tempLam-1)*tempLam^(2-j-z));           %replace this by a vector notation
        end
    end
%    toc(tic);       %testing
end

%% Start Tridiagonalization
indiag=zeros(bigL,1);
inrow=indiag;

inrow(2:bigL)  = gamma(1:bigL-1)./(2*sqrt(pi));
indiag(2:bigL) = xi(1:bigL-1);

[epsilon,t]=star2tridiag(indiag,inrow);

if modelpara.L == 0
    modelpara.L = find(epsilon > modelpara.precision, 1,'Last') + 1;    % +1 as L includes spin site
    dispif(sprintf('Optimum chain length is: %u',modelpara.L),modelpara.logging)
end

[epsilon,t]=extroplate(epsilon,t,modelpara.L);              % extrapolate very small levels for higher precision
modelpara.epsilon=epsilon(1:modelpara.L-1);
modelpara.t=t(1:modelpara.L-1);
%modelpara.t=abs(modelpara.t); %I noticed that the r form hess() function sometimes changed sign.
end