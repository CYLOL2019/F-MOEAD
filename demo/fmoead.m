function [pareto, apareto] = fmoead(mop)
%MOEAD run moea/d algorithms for the given mop.
    %global variable definition.
    global optimal_inds archive
    evolve(mop);
    pareto=[optimal_inds];
    apareto = [archive];
end

% The evoluation setp in MOEA/D
function evolve(mop)
    global params idealpoint nadirpoint cachepoints;
    selindex = 1:params.popsize;  
    selindex = random_shuffle(selindex);
    for i=1:params.popsize
        index = selindex(i);
        r = rand;
        useneighbour = r < params.neighbormating;
        ind = genetic_op_DEPM_swap(index,mop.domain,useneighbour);
        [ind] = getobjectives(ind);
        % tmp for ideal and nadir
        tmp = [cachepoints, ind.objectives];
        [idealpoint, I] = min(tmp,[],2);
        cachepoints = tmp(:,I);
        nadirpoint = max(cachepoints,[],2);
        update(index, ind, useneighbour);
    end
    % non-dominance and crowding distance, 1-by-1
%     archive = ndcd_1b1(optimal_inds,archive);
end
% update the population.
% index is the subproblem's index in the main population.
% ind is the individual structure.
% useneighbour is a bool determine whether the neighbourhood of index, or the whole population should be updated.
% this procedure is also governed by a parameter from params: params.updatesize, which determine how many subproblem
% should be updated at most by this new individual in MOEAD-DE.

%%
function update(index, ind, updateneighbour)
global idealpoint params optimal_inds;

% collect the updation index
if (updateneighbour)
updateindex = params.neighbor(1:params.Tr,index);
else
updateindex = [1:params.popsize]';
end

updateindex = random_shuffle(updateindex);
neighborpoints = [optimal_inds(updateindex)]; 


newobj  = subobjective([params.W(:,updateindex)], [ind.objectives], idealpoint, params.dmethod);    %objective values of current solution
oldobj  = subobjective([params.W(:,updateindex)], [neighborpoints.objectives], idealpoint, params.dmethod);        %previous objective values
C       = newobj < oldobj;    %new solution is better?
if (length(C(C==1)) <= params.updatesize)
    toupdate = updateindex(C);
else
    toupdate = randsample(updateindex(C), params.updatesize);
end  
% for i = toupdate'
[optimal_inds(toupdate')]=deal(ind);   %repace with the new one 
% end
end

%% get objectives 
function [x] = getobjectives(x)
% v = [];
for i  = 1 : length(x)
    vi   = cvar(x(i).combination);
    x(i).objectives = vi;

%     v = [v,vi];
end
end





