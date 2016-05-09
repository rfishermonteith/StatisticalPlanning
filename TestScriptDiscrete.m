%% Test script for discretised sampling

close all
clear
dispVis = 0;
n = 101;
r = 0.1;
tol = 1/n;

pStart = 0.65;
pEnd = 0.3;
probEnd = 0.05;

x = linspace(0,1,n);
p = ones(1,n)/n*(1-probEnd);
%p = normpdf(x, 0.5, 0.3)/n*(1-probEnd);
width = x(2)-x(1);
p(abs(x-pEnd)<tol) = p(abs(x-pEnd)<tol) + probEnd;

% overload val2ind:
val2ind = @(val, width) round(val*(n-1)+1);

tic

if dispVis
    figure;
    subplot(5,1,1)
    plot(x,p);
    ylim([0,1.1*max(p)])
    title('Sampling probabilities')
    
end

pNode{1} = zeros(1,n);
pNode{1}(x==pStart) = 1;

if dispVis
    subplot(5,1,2)
    plot(x,pNode{1});
    title('pNode 1')
end

% Run sim m times:
m = 200;
EX = 0;

% Precalculate upper and lower bounds:
upperBoundA = zeros(1, length(x));
lowerBoundA = zeros(1, length(x));
for kk = 1:length(x)
    lowerBoundA(kk) = val2ind(max(x(1),x(kk)-r), width);
    upperBoundA(kk) = val2ind(min(x(kk)+r,x(end)), width);
end

upperBoundB = zeros(length(x));
lowerBoundB = zeros(length(x));
for kk = 1:length(x)
    for kkk = 1:kk
        lowerBoundB(kk, kkk) = val2ind(max(x(1),x(kkk)-r), width);
        upperBoundB(kk, kkk) = val2ind(min(x(max(kk-1,1))+r,x(end)), width);
    end
end

upperBoundC = zeros(1, length(x));
lowerBoundC = zeros(1, length(x));
for kk = 1:length(x)
    for kkk = kk:length(x)
        lowerBoundC(kk, kkk) = val2ind(max(x(1),x(min(kk+1,length(x)))-r), width);
        upperBoundC(kk, kkk) = val2ind(min(x(kkk)+r,x(end)), width);
    end
end

for k = 2:m
    % Calculate or(p_(i-1)(x-r:x+r))
    orA = zeros(1, length(x));
    for kk = 1:length(x)
        orA(kk) = distOr2(pNode{k-1}(lowerBoundA(kk):upperBoundA(kk))');
    end
    
    % Calculate or(p(0:x)
    orB = zeros(1,length(x));
    for kk = 1:length(x)
        orB1 = zeros(1, length(x));
        for kkk = 1:kk-1 % DEBUG
            orB1(kkk) = (1-distOr2(pNode{k-1}(lowerBoundB(kk, kkk):upperBoundB(kk, kkk))'))*p(kkk);
        end
        orB(kk) = sum(orB1);
    end
    
    % Calculate or(p(x:end)
    orC = zeros(1,length(x));
    for kk = 1:length(x)
        orC1 = zeros(1,length(x));
        for kkk = kk+1:length(x) % DEBUG
            orC1(kkk) = (1-distOr2(pNode{k-1}(lowerBoundC(kk, kkk):upperBoundC(kk, kkk))'))*p(kkk);
        end
        orC(kk) = sum(orC1);
    end
    
    if dispVis
        subplot(5,1,3)
        hold off
        plot(x,orA);
        hold on
        ylim([0 1.1])
        plot(x,orB);
        plot(x,orC);
        legend('orA','orB','orC')
    end
    
    % Calculate shifts
    bound = find(abs(x - r)<tol);
    pLeft = [pNode{k-1}(bound:end), zeros(1,bound-1)];
    pRight = [zeros(1,bound-1), pNode{k-1}(1:end-bound+1)];
    
    if dispVis
        subplot(5,1,4)
        hold off
        plot(x, pLeft)
        hold on
        plot(x, pRight)
        legend('pLeft','pRight')
    end
    
    % Calculate tempNode{k}
    %tempProb = distOr2([p.*orA; pLeft.*orB; pRight.*orC]);
    tempProb = sum([p.*orA; pLeft.*orB; pRight.*orC]);
    
    %     for kk = 1:length(p)
    %         tempProb(kk) = distOr2([p(kk).*orA(kk); pLeft(kk).*orB(kk); (pRight(kk)*(1 - distOr2(pRight(kk+1:end)'))).*orC(kk)]);
    %     end
    % Normalise
    %tempProb = tempProb./sum(tempProb);
    pNode{k} = distOr2([tempProb;pNode{k-1}]);

    % Check that it sums to 1
    if abs(sum(tempProb)-1)>tol
        warning(['tempProb does not sum to 1, instead, it sums to ', num2str(sum(tempProb))])
    end
    if dispVis
        subplot(5,1,5)
        plot(x,tempProb);
    end
    
    if dispVis
        subplot(5,1,2);
        plot(x,pNode{k});
        title(['pNode ', num2str(k)])
        drawnow;
    end
    
    EX(k) = EX(k-1) + k*(pNode{k}(x==pEnd)-pNode{k-1}(x==pEnd));
    %     for k2 = 1:k-1
    %         probs(k2,:) = pNode{k2}(x==pEnd);
    %     end
    %     EX(k) = EX(k-1) + k*(pNode{k}(x==pEnd)*(1-prod(probs)));
    
    % Stop if prob at pEnd > 0.5
    %     if pNode{k}(abs(x-pEnd)<tol)>0.5
    %         disp(k)
    %         break
    %     end
    
end

figure; plot(EX);
toc

%% Run a simulation of standard RRT to test the results
tic

numIter = 10000;
numNodes = zeros(numIter,1);
for i = 1:numIter
    cont = true;
    count = 1;
    pointsList = pStart;
    
    while cont
        
        % Sample a point
        sample = generateRandomPoint(p,x, pEnd, 0.0);
        
        % Find nearest node
        distancesToSample = abs(pointsList - sample);
        [closestDistance, closestNode] = min(distancesToSample);
        direction = -(pointsList(closestNode) - sample)/abs(pointsList(closestNode) - sample);
        
        % If closestDistance is less than r, add it to pointsList, else move r toward it
        if closestDistance < r
            pointsList(end+1) = sample;
        else
            [~,indConnect] = min(abs(x-(pointsList(closestNode) + direction.*r)));
            pointsList(end+1) = x(indConnect);
            %pointsList(end+1) = round(pointsList(closestNode) + direction.*r,3);
        end
        
        if pointsList(end) == pEnd
            cont = false;
            numNodes(i) = length(pointsList);
        end
    end
end
expectedNodesRRT = mean(numNodes)
toc