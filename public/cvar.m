function pf = cvar(combination)
global mocpo_params
%% Parameters for AUGMECON2 CVaR
% parameters for the problem
roundlot = mocpo_params.Rl;
budget = 1/mocpo_params.Rl;
alpha = mocpo_params.alpha;
% vub = 1;
% vlb = 0.1;
ub = mocpo_params.Ub(combination)/mocpo_params.Rl;
lb = ceil(mocpo_params.Lb(combination)/mocpo_params.Rl);
grid_num = mocpo_params.grid;
obj_num = 2; % 1, CVaR 2, PR (Return)

% data
ret = mocpo_params.umx(combination,:)';
ti = mocpo_params.ti;
[intervals, noa] = size(ret);
ep = mean(ret); % 1,noa

% scenario probability
pr = ones(intervals,1)/intervals; 
pf = NaN(obj_num,grid_num);
ps = NaN(noa,grid_num);

% solver options
options = optimoptions('intlinprog','RelativeGapTolerance',0,'Display','off');
eps = 10 ^-5; % parameter in AUGMECON

% options = optimoptions('intlinprog','RelativeGapTolerance',0);
%% Construct Model
portprob = optimproblem;

% Variables
% portfolio = optimvar('portfolio',noa,'LowerBound',0);
portfolio = optimvar('portfolio',noa, 'Type', 'integer');
% sel = optimvar('sel',noa,'Type','integer','LowerBound',0,'UpperBound',1);
VaRdev = optimvar('VaRdev',intervals,'LowerBound',0);
VaR = optimvar('VaR');
losses = optimvar('losses',intervals);
CVaR = VaR + pr'*VaRdev/(1-alpha);
PR = - ep * portfolio * roundlot;

% Auxiliary Variables
totalbudget = sum(portfolio);

% Constraints
portprob.Constraints.portl = portfolio >= lb ; 
portprob.Constraints.portu = portfolio <= ub ; 
% portprob.Constraints.cardu = sum(sel) == Ku;
portprob.Constraints.conbudget = totalbudget == budget;
portprob.Constraints.conVaRdev = VaRdev >= losses - VaR;
portprob.Constraints.conlossDef = losses == (ti*budget - ret * portfolio)*roundlot;

% % Objective
% portprob.Objective = CVaR;

%% Get payoff table, reduce payoff bound in AUGMECON

% payoff = getminallobj(obj_num, portprob, CVaR, PR);
payoff = getminallobj(obj_num, portprob, CVaR, PR, ep, roundlot, pr, alpha);
pf(:,1) = payoff(1,:)';
pf(:,end) = payoff(2,:)';
slrange = payoff(1,2)-payoff(2,2);
% ps = NaN(noa,grid_num);

%% Iterative Optimization, jump duplicates in AUGMECON2
% if slrange == 0, it denotes the local PF is just one point
if slrange == 0
    pf(1,:) = pf(1,1);
    pf(2,:) = pf(2,1);
else
    %  slack or surplus variables for the eps-constraints
    sl = optimvar('sl',obj_num-1,'LowerBound',0);
    rhs  =  linspace(payoff(1,2),payoff(2,2), grid_num);
    step_gap = rhs(1) - rhs(2);
    % Objective
    portprob.Objective = CVaR + eps * sl/slrange;
    i = 2;
    while i < grid_num    
    %     i
        portprob.Constraints.conobjminor = PR == rhs(i) - sl;
        [sol,fval] = solve(portprob,'Options',options);   
        pf(1,i) = sol.VaR + pr'*sol.VaRdev/(1-alpha);
        pf(2,i) = rhs(i) - sol.sl;   
    %     ps(:,i) = sol.portfolio*roundlot;
        if sol.sl > 0 
            i = i + ceil(sol.sl/step_gap);
        else
            i = i + 1;
        end
    %     sol.sl/step_gap
    end
end




 
end


%% auxiliary functions
function v = getminallobj(obj_num, portprob, CVaR, PR, ep, roundlot, pr, alpha)
options = optimoptions('intlinprog','Display','off');
% options = optimoptions('intlinprog');
v = zeros(obj_num,obj_num);
for i = 1 : obj_num % rows
    for j = 1 : obj_num % columns
        switch i
            case 1
                switch j
                    case 1
                        portprob.Objective = CVaR;
                        portprob.Constraints.objcon = [];
                        [sol,fval] = solve(portprob,'Options',options);
                        v(1,1) = fval;
                        % variable tmp is used when intlinprog can not work
                        % well
                        tmp = - ep * sol.portfolio * roundlot;
                    case 2
                        portprob.Objective = PR;
                        portprob.Constraints.objcon = CVaR == v(1,1);
                        [sol,fval] = solve(portprob,'Options',options);
                        % some error with this intlinprog
                        if isempty(fval)
                            v(1,2) = tmp;
                        else
                            v(1,2) = fval;
                        end

                end
            case 2
                switch j
                    case 1
                        portprob.Objective = PR;
                        portprob.Constraints.objcon = [];
                        [sol,fval] = solve(portprob,'Options',options);
                        v(2,2) = fval;
                        % variable tmp is used when intlinprog can not work
                        % well
                        tmp = sol.VaR + pr'*sol.VaRdev/(1-alpha);
                    case 2
                        portprob.Objective = CVaR;
                        portprob.Constraints.objcon = PR == v(2,2);
                        [sol,fval] = solve(portprob,'Options',options);
                         % some error with this intlinprog
                        if isempty(fval)
                            v(2,1) = tmp;
                        else
                            v(2,1) = fval;
                        end
                end
        end
    end
end
end