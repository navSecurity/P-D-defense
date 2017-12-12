function [o] = simulatePincerParameters(s)
% simulatePincerParameters : Simulate the parameters in the vector theta, 
%                            on which observables depend, for
%                            hypothesis testing in the pincer defense.
%
%
% INPUTS
%
% s ------- structure variable with the following fields:
% 
%            N = Number of parameter sets to simulate
%
%            i = Hypothesis designator: 0 = null hypothesis, 1 = multipath, 2
%                = spoofing, 3 = jamming
%
%        elDeg = Elevation angle for the multipath scenario, in deg.  The
%                elevation angle can be thought of as a proxy for the severity
%                of multipath, with lower elevation being more severe.  
%
%           Tc = Chip interval (seconds)
%
%      muEtadB = Mean value of the power advantage marginal distribution for
%                multipath and spoofing (dB).  
%
%   sigmaEtadB = Std of the power advantage marginal distribution for
%                multipath and spoofing (dB). 
%
%      riceSdB = Rician s value of power advantage marginal distribution
%                for jamming (dB).
%
%  riceSigmadB = Rician sigma value of power advantage marginal
%                distribution for jamming (dB).
%
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
% Author:  Todd Humphreys
%+==============================================================================+

switch s.i
  case 0  
    o.etadB = -100;
    o.Delta_theta = zeros(s.N,1); 
    o.Delta_tau = zeros(s.N,1);    
  case 1
    s.maxEtadB = -1.5;
    % A minimum etadB is set to -40 to match the minimum signal strength of paths
    % in the LMSCM simulator.
    s.minEtadB = -40;  
    o = simulateMultipathParameters(s);
  case 2
    % Spoofing is assumed to mimic the multipath distribution but with higher, and
    % tighter, etadB, and with a broad range of delay, which can be simulated
    % as multipath from an SV at low elevation.
    s.minEtadB = -1.5;
    o = simulateMultipathParameters(s);
  case 3
    % Jamming will be assumed to have a Rician distribution for etadB.  This
    % becomes a Rayleigh distribution when the Rician parameter s = 0.
    pd = makedist('Rician','s',s.riceSdB,'sigma',s.riceSigmadB);
    o.etadB = random(pd,s.N,1);
    % The phase for jamming is arbitrary, as it has no effect on Bayes
    % risk.
    o.Delta_theta = 2*pi*rand(s.N,1);
    % The delay for jamming is arbitrary, so long as it exceeds two chips
    % so that there is no interference of correlation functions.
    o.Delta_tau = 3*ones(s.N,1);
   
end