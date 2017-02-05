function [o] = simulateMultipathParameters(s)
% simulateMultipathParameters : Simulate the parameters in the vector theta,
%                               on which observables depend, for hypothesis
%                               testing in the pincer defense under the
%                               multipath hypothesis.
%
%
% INPUTS
%
% s ------- structure variable with the following fields:
% 
%            N = Number of parameter sets to simulate
%
%        elDeg = Elevation angle for the multipath scenario, in deg.  The
%                elevation angle can be thought of as a proxy for the severity
%                of multipath, with lower elevation being more severe.  See
%                Figs. 12 to 19 in [1].  
%
%         view = Flag indicating whether to plot the joint and marginal
%                distributions of etaB and Delta_tau for inspection.
%
%           Tc = Chip interval (seconds)
%
%      muEtadB = Mean value of the power advantage marginal distribution (dB).
%                If this field is not present, the default value below is
%                used.
%
%   sigmaEtadB = Std of the power advantage marginal distribution (dB).  If
%                this field is not present, the default value below is
%                used.
%
%     maxEtadB = Maximum value of power advantage marginal distribution
%                (dB). Samples exceeding this value are replaced with other
%                samples that respect the bound. If this field is not present,
%                no maximum bound is enforced.
%
%     minEtadB = Minimum value of power advantage marginal distribution
%                (dB). Samples below this value are replaced with other
%                samples that respect the bound. If this field is not present,
%                no minimum bound is enforced.
%
% OUTPUTS
%
% o ------- structure variable with the following fields:
%
%        etadB = N-by-1 vector of MC samples of the power ratio between the
%                the most powerful multipath signal and the LOS signal (dB).
%
%    Delta_tau = N-by-1 vector of MC samples of tau_I - tau_A values, the
%                difference between the multipath and the authentic code
%                phases (chips)
%
%  Delta_theta = N-by-1 vector of MC samples of theta_I - theta_A values, the
%                difference between the multipath and the authentic carrier
%                phases (rad)
%
%+------------------------------------------------------------------------------+
% References: 
% 
% [1] Andreas Lehner et al., "Measuring the navigation multipath channel-a
%     statistical analysis," in Proceedings of the ION GNSS 2004.
%
% [2] Turin, George L., et al. "A statistical model of urban multipath
%     propagation." IEEE Transactions on Vehicular Technology 21.1 (1972):
%     1-9.
%
% [3] Andreas Lehner, "Multipath Channel Modelling for Satellite Navigation
%     Systems”, PhD Thesis, University of Erlangen-Nuremberg, ISBN
%     978-3-8322-6651-6, Shaker Verlag GmbH, D-52018 Aachen, Germany, 2007.
%
% Author:  Todd Humphreys
%+==============================================================================+

if(~isfield(s,'view'))
  s.view = 0;
end

% The default parameters for the statistical multipath model were derived from
% an analysis of the Land Mobile Multipath Channel Model described in [3].

% lognormal distribution for power advantage
if(isfield(s,'muEtadB'));
  muEta = s.muEtadB;
else
  muEta = -21;
end
if(isfield(s,'sigmaEtadB'));
  sigmaEta = s.sigmaEtadB;
else
  sigmaEta = 5;  
end
% exponential distribution for excess delay, with the exponential
% distribution parameter mu a quadratic function of elevation angle
polyCoeffs = [0.0120   -2.4061  133.7863]';
muDelay = polyval(polyCoeffs,s.elDeg);
% The linear correlation coefficient between the lognormal and exponential
% marginal distributions, rho, is approximately enforced by a three-step
% procedure: (1) generate a multivariate normal distribution with a linear
% correlation rho, (2) flatten the multivariate normal into a uniform
% distribution by passing through the normcdf, (3) shape the correlated
% uniform distribution into normal and exponential marginals with the
% correct parameters using norminv and expinv.
rhoEtaDelay = -0.26;
Z = mvnrnd([0 0], [1 rhoEtaDelay; rhoEtaDelay 1], 4*s.N);
U = normcdf(Z);
etadB = norminv(U(:,1),muEta,sigmaEta);
Delta_tau_ns = expinv(U(:,2),muDelay);
Delta_tau_ns = Delta_tau_ns(1:s.N);
% Enforce bounds on etadB
if(isfield(s,'maxEtadB'))
  maxEtadB = s.maxEtadB;
else
  maxEtadB = max(etadB);
end
if(isfield(s,'minEtadB'))
  minEtadB = s.minEtadB;
else
  minEtadB = min(etadB);
end
o.etadB = etadB(1:s.N);
spareEtadB = etadB(s.N + 1:end);
iidum = find((spareEtadB <= maxEtadB) & (spareEtadB >= minEtadB));
compliantEtadB = spareEtadB(iidum);
iidum = find((o.etadB > maxEtadB) | (o.etadB < minEtadB));
if(length(iidum) > length(spareEtadB))
  error('Not enough bound-compliant spare samples of eta');
end
o.etadB(iidum) = compliantEtadB(1:length(iidum));
% The carrier phase difference Delta_theta can be accurately modeled as
% independent from etadB and Delta_tau and uniformly distributed on
% [0,2*pi); see [2].
o.Delta_theta = 2*pi*rand(s.N,1);
% Convert Delta_tau_ns to chips
o.Delta_tau = 1e-9*Delta_tau_ns/s.Tc;

if(s.view)
  figure(1);clf;
  plot(o.Delta_tau_ns,o.etadB,'.');
  xlabel('\Delta \tau (ns)');
  ylabel('\eta (dB)');
  title('Scatter plot of power advantage \eta vs. delay \Delta\tau');
  figure(2);clf;
  hist(o.etadB,40);
  title('Histogram of \eta');
  xlabel('dB');
  figure(3);clf;
  hist(X(:,2),40);
  title('Histogram of \Delta\tau');
  xlabel('ns');
  fprintf('Linear correlation matrix:\n');
  R = corr([o.etadB,o.Delta_tau_ns])
end


