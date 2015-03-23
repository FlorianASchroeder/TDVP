function [modelpara,xi,Gamma]=SBM_genpara(modelpara)
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
%   FS 16/11/2014:  - Added analytical result for Orthogonal Polynomial
%                     chain mapping of J(w) of the SBM.
%	FS 12/02/2015:	- Adding Stieltje to OrthPol methods.
%					- Simplifying the conditionals for different methods
%					- Added more control over Analytic extraction and Numerical extraction of chain parameters.
%					  Allows for many more combinations!

s = modelpara.s;
% cutoffs are unlikely to be changed -> leave here!
if strcmp(modelpara.chain.spectralDensity, 'Renger')
	w_cutoff = cmToeV(3000);                            % upper cutoff frequency
else
	% All other cases use J(w) from Leggett -> SBM
	w_cutoff = 1;
end
% other declarations

% first: OrthPol, analytic from recurrence relation.
if strcmp(modelpara.chain.mapping,'OrthogonalPolynomials')
	assert(isfield(modelpara,'L'),'Please state required chain length!');
    assert(modelpara.L>1,'Please give finite value for chain length!');
    assert(isfield(modelpara,'alpha'),'Please state para.alpha, the strength of coupling!');

	if strcmp(modelpara.chain.method,'Analytic')
		[modelpara.epsilon, modelpara.t] = chainParams_OrthogonalPolynomials();
		return;																% Since everything is done!
	elseif strcmp(modelpara.chain.method,'Stieltjes')
		% Needs discretization, either Linear or Logarithmic!
		% Done in next stage!
	else
		error('Please define valid chain.method for Orthogonal Polynomials!');
	end
end

%% Start Discretization
% Used in combinations:
%	- OrthPol + Stieltjes
%	- LanzcosTriDiag
% either numerically using integrations (Renger)
% or from analytical result (Leggett)
% if para.Lambda > 1 -> Logarithmic, Zitko
% if para.Lambda = 1 -> Linear

z = modelpara.z;

if strcmp(modelpara.chain.discretization,'LogZ')
	assert(modelpara.Lambda ~= 1, 'Please define Lambda > 1 for use of LogZ discretization!');
	bigL = ceil(floor(-1*log(realmin)/log(modelpara.Lambda))/2) %Use a large enough start H to make sure the accuracy after transformed to Wilson chain
elseif strcmp(modelpara.chain.discretization,'Linear')
	assert(modelpara.L > 1,'You need to define a chain length > 1!');
	% TODO: make sensible estimate!
	bigL = 100*modelpara.L		% arbitrary now
else
	% 'None' is not applicable for following parts!
	error('You need to define a discretization method!');
end

if strcmp(modelpara.chain.discrMethod,'Analytic')
	assert(strcmp(modelpara.chain.discretization,'LogZ'), 'Analytic discretization only available for LogZ!');
	assert(strcmp(modelpara.chain.spectralDensity,'Leggett_Hard'), 'Analytic discretization only available for Leggett_Hard!');
	% If analytic discretization:
%% use J(w) from A. Leggett et al., Rev. Mod. Phys. 59, 1 (1987). DOI: 10.1103/RevModPhys.59.1
    %           and R. Bulla, N.-H. Tong, and M. Vojta, Phys. Rev. Lett. 91, 170601 (2003). DOI: 10.1103/PhysRevLett.91.170601
    %           and R. Žitko and T. Pruschke, Phys. Rev. B 79, 085106 (2009). DOI: 10.1103/PhysRevB.79.085106
    % Obtained by solving Eq. C4 in Žitko 2009, with f'(x) == 0;
    fac=2*pi*modelpara.alpha/(1+modelpara.s);
    tempfac=1/((1+modelpara.s)*log(modelpara.Lambda));
    tempLam=modelpara.Lambda^(1+modelpara.s);
    tempexp=1/(1+modelpara.s);

    xi=zeros(bigL,1);
    Gamma=zeros(bigL,1);

    % calculate discretized energy levels and couplings
	for j=2:bigL+1
		if j==2
            xi(j-1)=(tempfac*(1-tempLam^(-z))-z+1)^tempexp;
            Gamma(j-1)=sqrt(fac*(1-tempLam^(-z)));                      % modify this after the vector calculation.
        else
            xi(j-1)=(tempfac*tempLam^(2-j-z)*(tempLam-1))^tempexp;
            Gamma(j-1)=sqrt(fac*(tempLam-1)*tempLam^(2-j-z));           % replace this by a vector notation
		end
	end
else
	%% from here: either discretize J(w) for Lanzcos or h^2(x) for Stieltjes
	% Uses numerical evaluation of Integrals. Works for any spectral function!
	% numerical LogZ uses not the Zitko scheme, which is more difficult!!
% 	try
		% Following only works if spectral functions at the end have name scheme: J_(para.chain.spectralDensity)
		if strcmp(modelpara.chain.mapping,'LanzcosTriDiag')
% 			fprintf('J_%s(w,i)',modelpara.chain.spectralDensity)
			JorH = @(w,i) eval(sprintf('J_%s(w,i)',modelpara.chain.spectralDensity));
			wc = w_cutoff;
		elseif strcmp(modelpara.chain.mapping,'OrthogonalPolynomials')
			assert(strcmp(modelpara.chain.method,'Stieltjes'),'Discretization for Orthogonal Polynomials needs method Stieltjes!');
			JorH = @(x,i) eval(sprintf('J_%s(x.*w_cutoff,i)',modelpara.chain.spectralDensity))./pi.*w_cutoff;
			wc = w_cutoff;			% need to redefine cutoff = 1 since x = [0,1]!
		end
% 	catch
% 		error('cannot not find J_%s(w,i); Please define para.chain.spectralDensity properly!',modelpara.chain.spectralDensity);
% 	end

    xi	  = zeros(bigL,1);
    Gamma = zeros(bigL,1);			% not overloading gamma function!

	if strcmp(modelpara.chain.discretization,'LogZ')
		% gives different result than 'analytic', 'Leggett_Hard'
		% TODO: Problem: 0 is not included as boundary! need: w_limits = [w_limits,0];
		w_limits = modelpara.Lambda.^(1-z-(0:bigL)).*w_cutoff;      % Define limits of all the discretization intervals
		w_limits(1) = wc;
	elseif strcmp(modelpara.chain.discretization,'Linear')
		% Divide spectrum linearly
		% Problem: Soft cutoffs!
		w_limits = wc.*(1:-1/bigL:0);
	end

    for j=1:bigL
        intJ	  = integral(@(w) JorH(w,0),w_limits(j+1),w_limits(j));
        intJoverW = integral(@(w) JorH(w,1),w_limits(j+1),w_limits(j));
        xi(j)	  = intJ/intJoverW;							% bath energy levels
        Gamma(j)  = sqrt(intJ);								% bath couplings.
    end
end

xi(isnan(xi))=0;                                    % if numbers < e-161 they become NaN. This is fix for it.
Gamma(isnan(Gamma))=0;

%% Start mapping to Chain
if strcmp(modelpara.chain.mapping, 'LanzcosTriDiag')
	%% Start Tridiagonalization
	indiag=zeros(bigL,1);
	inrow=indiag;

	inrow(2:bigL)  = Gamma(1:bigL-1)./sqrt(pi);                         % factor of 1/2 only affected t(1) -> was moved into Hamiltonian; 1/sqrt(pi) changes only t(1)
	indiag(2:bigL) = xi(1:bigL-1);

	[epsilon,t]    = star2tridiag(indiag,inrow);

	if modelpara.L == 0
		modelpara.L = find(epsilon > modelpara.precision, 1,'Last') + 1;    % +1 as L includes spin site
		dispif(sprintf('Optimum chain length is: %u',modelpara.L),modelpara.logging)
	end

	[epsilon,t]			= extroplate(epsilon,t,modelpara.L);              % extrapolate very small levels for higher precision
	modelpara.epsilon	= epsilon(1:modelpara.L-1);
	modelpara.t			= t(1:modelpara.L-1);
	%modelpara.t=abs(modelpara.t); %I noticed that the r form hess() function sometimes changed sign.
elseif strcmp(modelpara.chain.mapping,'OrthogonalPolynomials') && strcmp(modelpara.chain.method,'Stieltjes')
	%% Start Stieltjes
	% assume: xi = levels = x
	%		  gamma = coupling strengths = h^2(x)
	ab = stieltjes(modelpara.L-1,[xi,Gamma.^2]);
	modelpara.epsilon = ab(:,1);
	modelpara.t		  = ab(:,2);
	% for eta_0: integrate [0,Inf]; any cutoff must be within the spectral function (stepfun)
	modelpara.t       = [sqrt(integral(@(w) eval(sprintf('J_%s(w,0)./pi',modelpara.chain.spectralDensity)),0,Inf));modelpara.t(1:end-1)];		% put sqrt(eta_0/pi) = t(1) by hand in!
else
	error('This is really hard to reach! (SBM_genpara)');
end
    function [w, t] = chainParams_OrthogonalPolynomials()
        %% returns w = epsilon and t for J(w) from A. Leggett
        % using orthogonal polynomial mapping from A. W. Chin et al. J. Math. Phys. 51, 092109 (2010). DOI: 10.1063/1.3490188
        % No need for tridiagonalization since analytical result.
        % t(1) == sqrt(eta_0/pi)/2
% 		wc = 1;     % cutoff is set to 1 always!
        n = 0:modelpara.L-2;
		if ~isempty(strfind(modelpara.chain.spectralDensity,'Leggett'))
			if ~isempty(strfind(modelpara.chain.spectralDensity,'Hard'))
				w = w_cutoff/2.*(1+ (s^2./((s+2.*n).*(2+s+2.*n))));                                     % w(1) = w(n=0)
				t = (w_cutoff.*(1+n).*(1+s+n))./(s+2+2.*n)./(3+s+2.*n).*sqrt((3+s+2.*n)./(1+s+2.*n));   % t(1) = t(n=0) != coupling to system
				t = [ w_cutoff*sqrt(2*pi*modelpara.alpha/(1+s))/(sqrt(pi)), t(1:end-1)];				% t(1) = sqrt(eta_0/pi)/2, 1/2 from sigma_z; t(2) = t(n=0)
			elseif ~isempty(strfind(modelpara.chain.spectralDensity,'Soft'))
				w = w_cutoff.*(2.*n+1+s);																% w(1) = w(n=0)
				t = w_cutoff.*sqrt((n+1).*(n+s+1));														% t(1) = t(n=0) != coupling to system
				t = [ w_cutoff*sqrt(2*pi*modelpara.alpha*gamma(s+1))/(sqrt(pi)), t(1:end-1)];           % t(1) = sqrt(eta_0/pi)/2, 1/2 from sigma_z; t(2) = t(n=0)
			end
		end
		w = w'; t = t';
	end

	function y = J_Renger(w,i)
		%% Spectral function for LH2 of Rhodoblastus acidophilus  from: F. Müh and T. Renger, Biochim. Biophys. Acta 1817, 1446 (2012). DOI: 10.1016/j.bbabio.2012.02.016
		% this code section can be used with any definition of J(w)
		% Units: eV
		% need for numerical integration.
		% multiplied by \pi * \omega^2: \pi to equal out division of \sqrt\pi later in the \gamma. Problem due to different
		%               definition of J(w) in Leggett et al.
		%               \omega^2 because Leggett defined J(w) as spectrum of the dimension-ful coupling constant! Renger
		%               took one \hbar\omega out.
		% This J(w) gives a Huang-Rhys factor of 1.3
		% i to get function J*w^(-i)

		%    tic; modelpara.Lambda = 2; bigL = 511; z=1;        % testing
		%   w_c = 1000; accurate up to 10^-3 eV; w_c = 2000: 10^-6 eV; w_c = 3000: 10^-8 eV
		% 	w_cutoff = cmToeV(3000);                            % upper cutoff frequency


% 		w=0:cmToeV(1):cmToeV(3000); i=0;                                                                % testing
		y = pi/(factorial(7)*2).*((0.8/cmToeV(0.56).^4 .* w.^(5-i).*exp(-(w/cmToeV(0.56)).^(1/2)))+(0.5/(cmToeV(1.94).^4) .* w.^(5-i).*exp(-(w/cmToeV(1.94)).^(1/2))));
% 		plot(y);xlabel('$\omega / cm^{-1}$');ylabel('$J(\omega)/eV$'); formatPlot(1);                             % testing

		% original definition: J = @(w) 1/(factorial(7)*2).*((0.8/cmToeV(0.56).^4 .* w.^3.*exp(-(w/cmToeV(0.56)).^(1/2)))+(0.5/cmToeV(1.94).^4 .* w.^3.*exp(-(w/cmToeV(1.94)).^(1/2))));
		% from Cheng:
		% J=@(w) 2*pi*modelpara.alpha*w.^modelpara.s;
	end

	function y = J_Leggett_Hard(w,i)
		%% Section for Spectral function for Spin-Boson Model from A. Leggett
		% Hard cutoff at wc using stepfunction
		% i to get function J*w^(-i)
% 		wc = 1;
		y = 2*pi*modelpara.alpha*w_cutoff.^(1-modelpara.s).*w.^(modelpara.s-i).*~stepfun(w,w_cutoff);
	end

	function y = J_Leggett_Soft(w,i)
		%% Section for Spectral function for Spin-Boson Model from A. Leggett
		% Soft cutoff at wc using exponential function
		% i to get function J*w^(-i)
% 		wc = 1;
		y = 2*pi*modelpara.alpha*w_cutoff^(1-modelpara.s).*w.^(modelpara.s-i).*exp(-w./w_cutoff);
	end

end