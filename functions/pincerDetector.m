function [ decision ] = pincerDetector( D, P, costType )
% pincerDetector: evaluate the P-D defense for a specified cost type
%
% INPUTS
%	D - symmetric difference magnitude 
%	P - power in dB
%	costType - 0: equally penalize for all misclassfication
% 	         - 1: fixed cost as listed in the paper
% 		 - 2: cost variable in theta as per the paper

if costType==2
load variableCostDecision.mat
elseif costType==1
 load fixedCostDecision.mat
elseif costType==0
  load fixedCostDecisionMisClass.mat
end
decision=interp2(deml,pow,pincerDecision,D,P,'nearest');

end

