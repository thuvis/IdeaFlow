function [ideaLabel1, ideaLabel2, ideaFlow_indicator, ideaFlow_leadLagTime] = ideaFlow(X, rank, threshold, clNum1, clNum2, timeClNum, taoMax)
% ideaFlow    Cluster words into ideas and aggregate word-level lead-lags into idea flows using a sparse tensor.
%
%   Input   : X                      : N1 x N2 x T x (2*taoMax +1) sparse tensor
%             rank                   : rank number of tensor factorization approximaton
%             threshold              : threshold for the proportion of correlated word pairs
%                                      between two ideas to decide whether they correlate
%             clNum1                 : word cluster or idea number of the first social group
%             clNum2                 : word cluster or idea number of the second social group
%             timeClNum              : time points cluster number
%             taoMax                 : maximum absolute value of lead-lag time (integer)
%
%   Output  : ideaLabel1             : N1 x 1 matrix. labels for each word of the first group
%                                      to tell which idea it is clustered into
%             ideaLabel2             : N2 x 1 matrix. labels for each word of the second group
%                                      to tell which idea it is clustered into
%             ideaFlow_indicator     : clNum1 x clNum2 x T matrix. iff (i, j, k) == 1, indicate the i-th idea of the first group
%                                      correlate with the j-th idea of the second group at the k-th time point
%             ideaFlow_leadLagTime   : clNum1 x clNum2 x T matrix. if (i, j, k) == t > 0, indicate the i-th idea of the first group
%                                      at the k-th time point flows to the j-th idea of the second group at the (k + |t|)-th time point;
%                                      if (i, j, k) == t <= 0, indicate the j-th idea of the second group at the k-th time point flows
%                                      to the i-th idea of the first group at the (k + |t|)-th time point;
%                                      if (i, j, k) == 0 in ideaFlow_indicator matrix, (i, j, k) in this matrix means nothing
%
%  Author   : Yangxin Zhong (zhongyx4869@163.com)

    warning off all;
    
    % Tensor factorization to cluster words into ideas.
    % Sometimes high rank will break down factorization algorithm. When this happens, lower down the rank.
    for r = rank:-1:1
        P = parafac_als(X, r); % Tensor factorization.
        if size(find(isnan(P.lambda)), 1) == 0
            break;
        end
    end
    
    % Use K-means to cluster the words of two social group, respectively.
    ideaLabel1 = kmeans(P.U{1}, clNum1);
    ideaLabel2 = kmeans(P.U{2}, clNum2);
    ideaFlow_indicator = zeros(clNum1, clNum2, size(X, 3));
    ideaFlow_leadLagTime = zeros(clNum1, clNum2, size(X, 3));
    
    % Tensor factorization to aggregate word-level lead-lags into idea flows.
    for i = 1:1:clNum1
        idea1 = find(ideaLabel1 == i);
        for j = 1:1:clNum2
            idea2 = find(ideaLabel2 == j);
            subTensor = X(idea1, idea2, :, :); % Extract subtensor limited to idea1 and idea2.
            if size(find(subTensor == 1), 1) == 0
                continue;
            end
            
            % Sometimes high rank will break down factorization algorithm. When this happens, lower down the rank.
            for r = rank:-1:1
                Q = parafac_als(subTensor, r);
                if size(find(isnan(Q.lambda)), 1) == 0
                    break;
                end
            end
            
            % Use K-means to cluster timepoints.
            timeCluster = kmeans(Q.U{end-1}, timeClNum);
            
            % Aggregate the consecutive timepoints in same cluster into lead-lag periods.
            timeRecord = 1;
            clusterRecord = timeCluster(1);
            for k = 2:1:(size(X, 3) + 1)
                if k == size(X, 3) + 1 || k ~= clusterRecord
                    % Extract subtensor limited to a certaion lead-lag period.
                    if ndims(subTensor) == 4
                        subTensor2 = subTensor(:, :, timeRecord:(k-1), :);
                        total = size(subTensor2, 1) * size(subTensor2, 2) * (k - timeRecord);
                    elseif ndims(subTensor) == 3
                        subTensor2 = subTensor(:, timeRecord:(k-1), :);
                        total = size(subTensor2, 1) * (k - timeRecord);
                    else
                        subTensor2 = subTensor(timeRecord:(k-1), :);
                        total = k - timeRecord;
                    end
                    
                    % Decide whether the two ideas have flow relationship in a certain period.
                    % If so, calculate the idea-level lead-lag time by averaging the word-level lead-lag times.
                    correlatedMatrix = find(subTensor2);
                    if size(correlatedMatrix, 1) >= threshold * total
                        leadLagTimeSet = correlatedMatrix(:, end);
                        leadLagTimeSet = leadLagTimeSet - (taoMax + 1); % Convert from row numbers to lead-lag times.
                        ideaFlow_indicator(i, j, timeRecord:(k-1)) = 1;
                        ideaFlow_leadLagTime(i, j, timeRecord:(k-1)) = mean(leadLagTimeSet);
                    end
                    
                    if k ~= size(X, 3) + 1
                        timeRecord = k;
                        clusterRecord = timeCluster(k);
                    end
                end
            end
        end
    end
    
end