function point = generateRandomPoint(p, x, pEnd, pEndPerc)

% pEndPerc of time sample pEnd
if rand < pEndPerc
    point = pEnd;
else


validSample = false;
while validSample == false;
% Generate random sample on x
choice = randi(length(x));
sample = x(choice);
%sample = (x(end)-x(1))*rand+x(1);
[~, closestPoint] = min(abs(x-sample));
if p(closestPoint) >= max(p)*rand
    validSample = true;
    point = x(closestPoint);
end
   
end
end