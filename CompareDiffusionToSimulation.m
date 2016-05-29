% Test whether the diffusion of TestScriptDiscrete2 matches the simulation

% This will be done by displaying the distribution of added nodes for the
% 1st N steps, as well as the cumulative node distribution
close all
clear

n = 101;
r = 0.1;
tol = 1/n;

pStart = 0.65;
pEnd = 0.0;
probEnd = 0.05;

x = linspace(0,1,n);
p = ones(1,n)/n*(1-probEnd);
%p = normpdf(x, 0.5, 0.3)/n*(1-probEnd);
%p = p + 1./(5-x);
%p = zeros(1, length(x));
%p(end) = 0.01;
width = x(2)-x(1);
p(abs(x-pEnd)<tol) = p(abs(x-pEnd)<tol) + probEnd;


% Normalise
p = p/sum(p);

% overload val2ind:
val2ind = @(val, width) round(val*(n-1)+1);

pNode{1} = zeros(1,n);
pNode{1}(x==pStart) = 1;

placedDist = p;

EX2 = 0;
EX3 = 0;
runningProduct = 1;
runningProduct3 = 1;

simEX = 0;

figure;

N = 201;
M = 100000;
simNodes = zeros(M,N);
simNodes(:,1) = repmat(pStart, M, 1);

allNodeDists = zeros(N,length(x));
allNodeDists(1,:) = pNode{1};

% simNodes = mat2cell(repmat(pStart,M, 1), ones(M,1));

for k = 2:N
    %% Do RRT simulation for M examples
    for kk = 1:M
        % Get current pointsList
        pointsList = simNodes(kk,:);
        pointsList = pointsList(1:k-1);
        % Sample a point
        sample = generateRandomPoint(p,x, pEnd, 0.0);
        
        % Find nearest node
        distancesToSample = abs(pointsList - sample);
        [closestDistance, closestNode] = min(distancesToSample);
        direction = -(pointsList(closestNode) - sample)/abs(pointsList(closestNode) - sample);
        
        % If closestDistance is less than r, add it to pointsList, else move r toward it
        if closestDistance <= r
            pointsList(k) = sample;
        else
            [~,indConnect] = min(abs(x-(pointsList(closestNode) + direction.*r)));
            pointsList(k) = x(indConnect);
        end
        
        simNodes(kk,1:k) = pointsList;
        
    end
    
    % Find the liklihood of having solved the problem on this iteration:
    for i = 1:M
        solved(i) = find(simNodes(i,:) == pEnd,1);
    end
    sumSolved = sum(solved == k)./M;
    simEX(k) = simEX(k-1) + k*sumSolved;
    
    % Plot the placed node positions
    subplot(2,1,1)
    countsPlaced = zeros(1, length(x));
    for kk = 1:length(x)
        countsPlaced(kk) = sum(simNodes(:,k)==x(kk));
    end
    plot(x, countsPlaced./sum(countsPlaced), 'b')
    
    % Plot the new cumulative of the nodes
    subplot(2,1,2)
    countsCumulative = zeros(1, length(x));
    for kk = 1:length(x)
        countsCumulative(kk) = sum(sum(simNodes(:,1:k)==x(kk), 1));
    end
    plot(x, countsCumulative, 'b')
    drawnow;
    
    %% Run the step of the forward method
    
    movedDist = zeros(1, length(x));

    % Loop through each index, and move it:
    for kk = 1:length(x)
        % Where should p[kk] move to?
        
        % Is there a node within distance r?
        lowerBound = val2ind(x(kk)-r, width);
        if lowerBound < 1, lowerBound = 1; end
        upperBound = val2ind(x(kk)+r, width);
        if upperBound > length(x), upperBound = length(x); end
        pWithinR = distOr2(pNode{k-1}(lowerBound:upperBound));
        
        movedDist(kk) = movedDist(kk) + pWithinR*placedDist(kk);
        
        % Where is the nearest node (outside of r)?
        pOfBeingClosestNode = zeros(1, length(x)); % probability of this being the nearest node to kk
        %start at kk, and move outwards, updating the probabilities as you
        %go
        pOfBeingClosestNode(kk) = pNode{k-1}(kk);
        pNodeFound = pOfBeingClosestNode(kk);
        for kkk = 1:length(x) % length(x) is the furthest we could go
            % What is the probability of a node at kk-kkk?
            lowerBound = kk-kkk;
            if lowerBound < 1, pLower = 0; else pLower = pNode{k-1}(lowerBound); end
            upperBound = kk+kkk;
            if upperBound > length(x), pUpper = 0; else pUpper = pNode{k-1}(upperBound); end
            
%             % DEBUG: scale pNode for pLower, so that the remaining
%             % possible points sum to 1
%             indicesToExclude = (kk-kkk:kk+kkk);
%             tempSums = sum(allNodeDists(1:k-1,[1:kk-kkk,kk+kkk:end]),2);
%             tempSums(tempSums == 0) = 1; % When this happens, the node from this row has already been used, so the only thigs which will be scaled will be zeros
%             pTempDist = bsxfun(@rdivide, allNodeDists(1:k-1,:),tempSums);
%             if size(pTempDist,1) ~= 1
%                 pTempDist = distOr2(pTempDist);
%             end
%             if pLower ~= 0
%                 pLower = pTempDist(lowerBound);
%             end
%             if pUpper ~= 0
%                 pUpper = pTempDist(upperBound);
%             end
            
            chosenLower = 0;
            if pLower ~= 0 && lowerBound >= 1
                chosenLower = pLower;
                if pUpper ~=0
                    % Probability lower but not upper:
                    p1 = pLower*(1-pUpper);
                    % Probability neither lower nor upper:
                    p2 = (1-pLower)*(1-pUpper);
                    % Probability upper but not lower:
                    p3 = pUpper*(1-pLower);
                    % Probability pLower and pUpper
                    p4 = pLower*pUpper;
                    
                    pCumulative = p1 + (pLower/(pUpper+pLower))*p4;
                    pOfBeingClosestNode(lowerBound) = pCumulative*(1-pNodeFound);
                    chosenLower = pLower;
                    %                     pOfBeingClosestNode(lowerBound) = (pLower/(pUpper+pLower))*(1-pNodeFound);
                else
                    pOfBeingClosestNode(lowerBound) = pLower*(1-pNodeFound);
                end
            end
            
            chosenUpper = 0;
            if pUpper ~= 0 && upperBound <= length(x) % DEBUG
                chosenUpper = pUpper;
                if pLower ~= 0
                    % Probability lower but not upper:
                    p1 = pLower*(1-pUpper);
                    % Probability neither lower nor upper:
                    p2 = (1-pLower)*(1-pUpper);
                    % Probability upper but not lower:
                    p3 = pUpper*(1-pLower);
                    % Probability pLower and pUpper
                    p4 = pLower*pUpper;
                    
                    pCumulative = p3 + (pUpper/(pUpper+pLower))*p4;
                    pOfBeingClosestNode(upperBound) = pCumulative*(1-pNodeFound);
                    chosenUpper = pUpper;
                else
                    pOfBeingClosestNode(upperBound) = pUpper*(1-pNodeFound);
                end
            end
            
            pNodeFound = distOr2([pNodeFound, pUpper, pLower]);
            
        end
        if abs(sum(pOfBeingClosestNode)-1) > 1e-6
            warning(['pOfBeingClosestNode sums to ',num2str(sum(pOfBeingClosestNode)), ', rather than 1'])
        end
        
        
        % Move this node toward pOfBeingClosestNode (shifted by r)
        indR = val2ind(r, width);
        
        for kkk = 1:length(x)
            if kkk <= kk-indR % Move node left
                movedDist(kkk+indR-1) = movedDist(kkk+indR-1) + placedDist(kk)*pOfBeingClosestNode(kkk);
            elseif kkk >= kk+indR % Move node right
                movedDist(kkk-indR+1) = movedDist(kkk-indR+1) + placedDist(kk)*pOfBeingClosestNode(kkk);
            else % Place node where it is (with whatever frequency there's a not within distance r
                % Actually do nothing for now, this should be handle above
            end
            
            
        end
        % Ensure that all of p(kk) has been placed
        if abs(sum(movedDist) - sum(placedDist(1:kk))) > 10e-6
            warning(['Not all of placedDist(',num2str(kk),') has been placed'])
        end
    end
    
    pNode{k} = distOr2([movedDist;pNode{k-1}]);
    pNode{k}(x==pEnd) = 0; % DEBUG!
    allNodeDists(k,:) = movedDist;
    % scale: DEBUG
    %scaleFactorThing = (1-length(x))/(sum(pNode{k})-length(x))
    %pNode{k} = 1-(1-pNode{k}).*scaleFactorThing;
   
    
    EX2(k) = EX2(k-1) + k*(movedDist(x==pEnd))*runningProduct;
    %EX3(k) = EX3(k-1) + k*(counts{k-1}(x==pEnd))./sum(counts{k-1})*runningProduct3;
    %runningProduct3 = runningProduct3*(1-(counts{k-1}(x==pEnd))./sum(counts{k-1}));
    runningProduct = runningProduct*(1-(movedDist(x==pEnd)));
    
    % Plot the distribution of new nodes
    subplot(2,1,1)
    hold on
    plot(x, movedDist,'r')
    hold off
    mse = mean((movedDist-countsPlaced./sum(countsPlaced)).^2);
    title(['Probability of placing a node, with mse: ',num2str(mse)])
    
    subplot(2,1,2)
    hold on
    plot(x, sum(allNodeDists)*M, 'r')
    hold off
    drawnow;
    % Plot the new cumulative node positions
end
