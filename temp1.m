numNodes = zeros(10,1);
for i = 1:10
    cont = true;
    count = 1;
    pointsList = pStart;
    
    while cont
        
        % Sample a point
        sample = generateRandomPoint(p,x);
        
        % Find nearest node
        distancesToSample = abs(pointsList - sample);
        [closestDistance, closestNode] = min(distancesToSample);
        direction = -(pointsList(closestNode) - sample)/abs(pointsList(closestNode) - sample);
        
        % If closestDistance is less than r, add it to pointsList, else move r toward it
        if closestDistance < r
            pointsList(end+1) = sample;
        else
            pointsList(end+1) = pointsList(closestNode) + direction.*r;
        end
        %     figure;
        %     histogram(pointsList);
        %     title(['RRT after ', num2str(length(pointsList)), ' iterations'])
        
        if pointsList(end) == pEnd
            cont = false;
            numNodes(i) = length(pointsList)
        end
    end
end
expectedNodesRRT = mean(numNodes)