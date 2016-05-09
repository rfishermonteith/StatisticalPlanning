function newDist = distOr2(probabilities)
% DISTOR calculates p(p1 U p2 U ...) for a discrete distribution
if length(probabilities) == 1
    newDist = probabilities;
elseif isempty(probabilities)
    newDist = 0;
else
    newDist = 1 - prod(1-probabilities);
end




