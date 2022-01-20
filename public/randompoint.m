function ind = randompoint(prob, n)
%RANDOMNEW to generate n new point randomly from the mop problem given.
global mocpo_params
if (nargin==1)
    n=1;
end
% randarray = rand(prob.pd, n);
randarray = rand(prob.pd, n);
lowend = prob.domain(:,1);
span = prob.domain(:,2)-lowend;
point = randarray.*(span(:,ones(1, n)))+ lowend(:,ones(1,n));
cellpoints = num2cell(point, 1);
% combinations for MP_gams
selpoint = decode(point,mocpo_params.PreAss,mocpo_params.K);
selcellpoints = num2cell(selpoint, 1);
indiv=get_structure('individual_MOCPO');
ind = repmat(indiv, 1, n);
[ind.parameter] = cellpoints{:};
[ind.combination] = selcellpoints{:};
end