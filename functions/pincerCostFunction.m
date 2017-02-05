function [ cost ] = pincerCostFunction( decision,...
                                        trueScenario,...
                                        etadB,...
                                        Delta_tau,...
                                        Delta_theta,...
                                        tau_eml,...
                                        costType)
% pincerCostFunction: return cost of deciding 'decision' given that
%                     'trueScenario' is the true simulation scenario.
%
% INPUTS
%
%         decision = declared detection decision
%                    1:Null, 2:Multipath, 3:Spoofing, 4:Jamming
%
%     trueScenario = the true simulated scenario
%                    1:Null, 2:Multipath, 3:Spoofing, 4:Jamming
%
% %%%%%%%  Elements of theta %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%            etadB = Power advantage of interference signal (dB)
%
%        Delta_tau = tau_I - tau_A, the difference between the interference
%                    and the authentic code phases (chips)
%
%      Delta_theta = theta_I - theta_A, the difference between the
%                    interference and the authentic carrier phase (rad)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%          tau_eml = Early-minus late correlator spacing used for code phase
%                    tracking, in chips
%
%         costType = 0 (default): uniformly penalize misclassification
%                                 C=ones(4,4)-eye(4)
%                    1: apply a fixed cost matrix as defined in code below
%                    2: apply a theta-dependent cost
%
% OUTPUTS
%             cost = cost of choosing the input decision given, normalized
%                    to reside on [0,1].
%
%+------------------------------------------------------------------------------+
% References:
%
% Author:  Jason Gross and Todd Humphreys
%+==============================================================================+

if(nargin < 7)
    costType = 0;
end

% Note:  Matlab indexing used throughout

% Explanation of cost matrix for costType = 1 (fixed-cost matrix)
%
% C[i,i] = 0 for all i:  correct decision
%
% C[1,2] = 0.2: low cost: some multipath-induced code- and carrier- phase
%               error might have been mitigated via receiver's multipath
%               mitigation routines had multipath been detected.
%
% C[1,3] = 1.0: high cost: spoofing goes undetected; receiver potentially
%               reports hazardously misleading information to user or
%               downstream subsystems.
%
% C[1,4] = 0.9: high cost: jamming goes undetected; receiver potentially
%               reports hazardously misleading information to user or
%               downstream subsystems.  Not so costly as C[1,3] because a
%               jammer has less control over receiver output than a
%               spoofer.
%
% C[2,1] = 0.1: low cost: receiver wastes computational resources trying to
%               mitigate phantom multipath. Receiver also runs a chance of
%               biasing code and carrier measurements when attempting to
%               excise phantom multipath contributor.
%
% C[2,3] = 1.0: high cost: spoofing goes undetected; receiver potentially
%               reports hazardously misleading information to user or
%               downstream subsystems.  Multipath mitigation deployed by
%               receiver can't be assumed to reduce this cost.
%
% C[2,4] = 0.9: high cost: jamming goes undetected; receiver potentially
%               reports hazardously misleading information to user or
%               downstream subsystems.  Not so costly as C[2,3] because a
%               jammer has less control over receiver output than a
%               spoofer. Multipath mitigation deployed by receiver can't be
%               assumed to reduce this cost.
%
% C[i,j] = 0.4 for all i in {3,4} and j in {1,2}: moderate cost: the cost of a
%              false alarm and of broken continuity in navigation solution
%              availability.
%
% C[4,3] = C[3,4] = 0.2: low cost: although the receiver misclassified
%                        spoofing as jamming or vice-versa, it correctly
%                        prevented hazardously misleading information by
%                        making its navigation solution unavailable
% 
CF = [0.0  0.2  1.0  0.9;
      0.1  0.0  1.0  0.9;
      0.4  0.4  0.0  0.2;
      0.4  0.4  0.2  0.0];

% Code phase error saturation value: full cost for code phase errors in excess
% of this value.  Obviously, this value is application dependent.
tauErrorSaturationChips = 0.3;
% Saturation cost for severe multipath
severeMultipathSaturationCost = 0.8;
% Saturation value of eta under jamming: full cost for eta above this value.
jamming_etaSaturationdB = 10;

if costType == 0
  % Uniform misclassification cost
  C = ones(4,4) - eye(4);
elseif costType == 1
  % Fixed cost matrix 
  C = CF;
else 
  % theta-dependent cost begins with CF as baseline
  C = CF;
  
  % C[1,2] : Cost depends linearly on the severity of error in code phase
  %          induced by multipath up to tauErrorSaturationChips, at which
  %          point the cost tops out at severeMultipathSaturationCost.  If
  %          necessary, the formula below can be modified to penalize carrier
  %          phase errors for applications such as CDGNSS that depend on
  %          carrier phase.
  if(decision == 1 && trueScenario == 2)
    [tauError,thetaError] = observationError(etadB,Delta_tau,Delta_theta);
    C(1,2) = min(severeMultipathSaturationCost,tauError/tauErrorSaturationChips);
  end
  
  % C[1,3] : Full cost in CF[1,3] if spoofing attack could capture the target
  %          receiver; otherwise, treat just as C[1,2]: the spoofing attack's
  %          effect is no worse than multipath.
  if(decision == 1 && trueScenario == 3)
    if(Delta_tau < tau_eml/2 && etadB < -1)
      [tauError,thetaError] = observationError(etadB,Delta_tau,Delta_theta);
      C(1,3) =  min(severeMultipathSaturationCost,tauError/tauErrorSaturationChips);
    end
  end
  
  % C[1,4], C[2,4] : Cost depends linearly on etadB up to
  %                  jamming_etaSaturationdB, beyond which full cost in
  %                  CF[1,4] is applied.
  C(1,4) = min(CF(1,4),etadB/jamming_etaSaturationdB);
  C(2,4) = C(1,4);
   
  % C[2,3] : Full cost in CF[2,3] if spoofer could capture receiver's
  %          tracking loops; otherwise, zero cost:  the spoofing attack's
  %          effect is no worse than multipath and the scenario was
  %          classified as multipath.  
  if(Delta_tau < tau_eml/2 && etadB < -1)
    C(2,3) = 0;
  end
  
end

cost = C(decision,trueScenario);






