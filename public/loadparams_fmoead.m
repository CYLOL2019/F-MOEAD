function params = loadparams_fmoead()
%LOADPARMS_ICDV init the parameters for F-MOEA/D
params = get_structure('fmoead_parameter');
% General Setting for MOEA/D
%params.seed = 177; % random seed.
% rng(0);
params.dmethod = 'te'; % decomposition method of choice.
params.F = 0.5; % the F rate in DE. %%%%params.F = [0.8 0.6 0.5 ];
params.CR = 0.9; % the CR rate in DE. %%%%params.CR = [0.2 0.9 1.0 ];
params.neighbormating = 0.9; % percentage to update the neighbour or the whole population
params.updatesize = 2; % the maximum updation number.
params.popsize = 100; % popsize
% params.evaluation = 10^5; % maximum evaluations
% params.offsprings = 10^5/ceil(sqrt(params.popsize))/2; % maximum offsprings
params.runtime = 2*10^3;
end

