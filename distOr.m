function newDist = distOr(p1, p2, width)
% DISTOR calculates p(p1 U p2) for a discrete distribution
% DEBUG:
width = 1;
newDist = p1 + p2 - p1.*p2.*width;


