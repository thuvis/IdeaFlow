%  A simple example to show how to use function ideaFlow.
%  Author   : Yangxin Zhong (zhongyx4869@163.com)

M = load('test.txt'); % Load the index matrix of the sparse tensor.
X = sptensor(M(2:end, 1:end),1,[M(1,1) M(1,2) M(1,3) M(1,4)]); % Convert the index matrix into a sparse tensor.

% Use ideaFlow function to discover ideas and their flow relationshps.
[ideaLabel1, ideaLabel2, ideaFlow_indicator, ideaFlow_leadLagTime] = ideaFlow(X, 4, 0.6, 2, 3, 3, (M(1,4)-1)/2);

fprintf('\n');

% Set the following three virables to check the result of different word pairs (i, j) at different timepoints (k).
i = 1; j = 2; k = 3;
% Now you can check which idea a certain word belongs to by ideaLabel1 (for words in social group 1) and ideaLabel2 (for words in social group 2).
disp(['Word', num2str(i), ' in social group1 belongs to idea', num2str(ideaLabel1(i)), ' of group1.']);
disp(['Word', num2str(j), ' in social group2 belongs to idea', num2str(ideaLabel2(j)),  ' of group2.']);
fprintf('Idea%d of group1 consists of word:', ideaLabel1(i));
fprintf('%3d', find(ideaLabel1 == ideaLabel1(i)));
fprintf('.\n');
fprintf('Idea%d of group2 consists of word:', ideaLabel2(j));
fprintf('%3d', find(ideaLabel2 == ideaLabel2(j)));
fprintf('.\n');

% And you can also check whether two ideas have a flow relationship at a certain timepoint by ideaFlow_indicator.
if ideaFlow_indicator(ideaLabel1(i), ideaLabel2(j), k) == 1
    disp(['Idea', num2str(ideaLabel1(i)), ' of group1 has flow relationship with idea', num2str(ideaLabel2(j)),  ' of group2 at timepoint ', num2str(k), '.']);
    % If they do have flow relationship, you can check the idea-level lead-lag time between them at the timepoint by ideaFlow_leadLagTime.
    llt = ideaFlow_leadLagTime(ideaLabel1(i), ideaLabel2(j), k);
    if llt > 0
        disp(['Idea', num2str(ideaLabel1(i)), ' of group1 flows to idea', num2str(ideaLabel2(j)),  ' of group2 with a lead time ',  num2str(llt), ' at this period.'])
    else
        disp(['Idea', num2str(ideaLabel2(j)), ' of group2 flows to idea', num2str(ideaLabel1(i)),  ' of group1 with a lead time ',  num2str(-llt), ' at this period.'])
    end
else
    disp(['Idea', num2str(ideaLabel1(i)), ' of group 1 has no flow relationship with idea', num2str(ideaLabel2(j)),  ' of group2 at timepoint ', num2str(k), '.']);
end