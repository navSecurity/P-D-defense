function decision  = pdmlDetector(D, P, costType)
% pdmlDetector: evaluate the PD-ML defense for a specified cost type
%
% INPUTS
%	D - symmetric difference magnitude 
%	P - power in dB
%	costType - 0: equally penalize for all misclassfication
% 	         - 1: fixed cost as listed in the paper
% 		 - 2: cost variable in theta as per the paper
% Author: Jason N. Gross


if costType==2
    load detectionCF2a1.mat
elseif costType==1
    load detectionCF1a1.mat
elseif costType==0
    load detectionCF0a1.mat
end
decision=interp2(d1,pow,pdmlDecision,D,P,'nearest');

end

