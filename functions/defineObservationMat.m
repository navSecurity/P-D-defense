function [ H ] = defineObservationMat( delays, taps, tapsize)
%% define observation matrix E14 from Blanco-Delgado
% Author: Jason N. Gross
H=zeros(length(taps),length(delays));
for j=1:length(delays)
    for i = 1:length(taps)
        H(i,j)=GPSCAAutoCorr((taps(i)+delays(j))*tapsize);
    end
end

