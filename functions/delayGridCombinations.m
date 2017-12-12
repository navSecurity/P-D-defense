function [ combos ] = delayGridCombinations(delays )
% only for upo to 2 rays for now
% Author: Jason N. Gross
offset=0;
k=0;
for t0 = 1:length(delays)
    for t1 = 1:length(delays)
        if(delays(t0) < (delays(t1)-offset))
            k=k+1;
            combos(k,:)=[delays(t0) delays(t1)]';
        end
    end
    
end

