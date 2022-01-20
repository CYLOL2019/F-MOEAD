%%
function [archive] = ndcd_1b1(population, maxPop)


%% NSGA-II Parameters
[nObj,nPop] = size(population.objectives);        % Population Size
% nVard = size(population.parameter,1);
% parameter = zeros(nVard,maxPop);
% objective = zeros(nObj,maxPop);
%% Initialization
empty_individual.parameter=[];
empty_individual.weights=[];
%
% empty_individual.Position=[];
empty_individual.Cost=[];
empty_individual.Rank=[];
empty_individual.DominationSet=[];
empty_individual.DominatedCount=[];
empty_individual.CrowdingDistance=[];


pop    = repmat(empty_individual,nPop,1);

for i=1:nPop
    
    pop(i).parameter = population.parameter(:,i);
    pop(i).Cost = population.objectives(:,i);

end

indiv=get_structure('individual_MultiScale');

% Merging
pop=[pop];  

% Non-Dominated Sorting
[pop, F]=NonDominatedSorting(pop);

if length(F{1})<=maxPop
    archive = repmat(indiv, 1, length(F{1}));
    for i = 1 : length(F{1})
        archive(i).parameter = pop(i).parameter;
        archive(i).weights = pop(i).weights;
        archive(i).objectives = pop(i).Cost;    
    end
else
    
    archive = repmat(indiv, 1, maxPop);
   
    
    pop_eliminate = pop(F{1});
    while (length(pop_eliminate) > maxPop)
        F_eliminate = 1:length(pop_eliminate);
        pop_eliminate = CalcCrowdingDistance(pop_eliminate,{F_eliminate});
        % Sort Based on Crowding Distance
        [~, CDSO]=sort([pop_eliminate.CrowdingDistance],'descend');
        pop_eliminate(CDSO(end))=[];
    end
    for i = 1 : maxPop
        archive(i).parameter = pop_eliminate(i).parameter;
        archive(i).weights = pop_eliminate(i).weights;
        archive(i).objectives = pop_eliminate(i).Cost;    
    end
end
    







end



