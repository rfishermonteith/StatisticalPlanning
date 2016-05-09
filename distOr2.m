function newDist = distOr2(probabilities)
% DISTOR calculates p(p1 U p2) for a discrete distribution
if length(probabilities) == 1
    newDist = probabilities;
elseif isempty(probabilities)
    newDist = 0;
else
    newDist = probabilities(1,:);
    for k = 2:size(probabilities,1)
        newDist = newDist + probabilities(k,:) - newDist.*probabilities(k,:);
    end
end




