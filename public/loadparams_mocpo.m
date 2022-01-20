function [mocpo_params] = loadparams_mocpo(testname,dirtime,pn)
global start_marketime end_marketime params
%LOADPARMS init the parameters for F-MOEA/D
csvdir = sprintf('../%s',dirtime);
% marketime = '20201130_train';
% Parameters for MOCPO (constraints)
% mocpo_params.Kl = 4;
mocpo_params.K = 10;
mocpo_params.Rl = 0.008;
mocpo_params.PreAss = [30];
mocpo_params.alpha = 0.99; % \alpha = 99% for CVaR
mocpo_params.grid = ceil(sqrt(params.popsize)); % grids for AUGMECON
% Since 1 is divisible for Rl, it can be regarded as mix integer
% programming
% load data from csv
csvfile = sprintf('%s/port%d_umx.csv',csvdir,pn);
tmp = readtable(csvfile, 'HeaderLines',1); 
mocpo_params.umx = tmp(:,2:end);
mocpo_params.c_VaRDev = size(mocpo_params.umx,1) * (1-mocpo_params.alpha);
mocpo_params.alphaT = ceil(mocpo_params.c_VaRDev);
mocpo_params.umx = mocpo_params.umx{:,:};
mocpo_params.umx = mocpo_params.umx';%Asset, Time
mocpo_params.umx = mocpo_params.umx + 1;

mocpo_params.NoA = size(mocpo_params.umx,1);
mocpo_params.u = mean(mocpo_params.umx,2);

% TargetIndex       //= 1.01 * ones(intervals,1);
csvfile = sprintf('%s/000300.SH_20201211.csv',csvdir);
T = readtable(csvfile,'PreserveVariableNames', true);
ti = (T{:,2:3});
date_index =intersect(find(ti(:,1)>=str2num(start_marketime)),find(ti(:,1)<=str2num(end_marketime)));
ti = ti(date_index,2);
mocpo_params.ti = 1 + (ti - ti(end))/ti(end);

mocpo_params.Lb = 0.01*ones(mocpo_params.NoA,1);
mocpo_params.Ub = ones(mocpo_params.NoA,1);
% data for swap operators
%Highest return
[~, mocpo_params.u_idx] = sort(mocpo_params.u, 'descend');
%Lowest Epexceted loss of t time invervals (related to CVaR)
[tmp, ~]     = sort(mocpo_params.umx-mocpo_params.ti', 2,'ascend');%low return part
[~, mocpo_params.cvar_idx] = sort(mean(tmp(:,1:mocpo_params.alphaT),2),'ascend');
%Lowest loss at time T (related to VaR)
[~, mocpo_params.var_idx] = sort(tmp(:,mocpo_params.alphaT),'ascend');


% % Normalization for the AUGMECON2&CPLEX
% norfile = sprintf('%s/payoff_%d,%d.mat',csvdir,start_marketime,end_marketime);
% tmp = load(norfile);%risk -return, min risk and max return
% tmp = tmp.payoff;
% tmp(:,2) = -tmp(:,2);
% [mu,~] = sort(mocpo_params.u,'descend');
% % Least multipliers that y*Rl>=Lb
% least_multipliers = ceil(mocpo_params.Lb/mocpo_params.Rl);
% maxweight = mocpo_params.SumOne-least_multipliers*(mocpo_params.K-1);
% weight = [maxweight; least_multipliers*ones(mocpo_params.K-1,1)]*mocpo_params.Rl;
% tmp(2,2) = mu(1:mocpo_params.K)'*weight;
% % normalization scalar
% portnor=[tmp(2,1),tmp(1,1);tmp(2,2),tmp(1,2)];
% [NULL,mocpo_params.objref] = mapminmax(portnor,0,1);
% mocpo_params.nor(3) = portnor(2,1);
% mocpo_params.nor(4) = portnor(2,2);
% mocpo_params.nor(1) = portnor(1,1);
% mocpo_params.nor(2) = portnor(1,2);
% testmop
testmop(testname, mocpo_params.NoA);% MOP could be obtained by function 'testmop.m'.

%   % gamsparams wgdx
% gamsparams.i.name = 'i';
% gamsparams.i.uels = [];
% for jm = 1 : mocpo_params.K
%     uelstr = sprintf('%d',jm);
%     gamsparams.i.uels = [gamsparams.i.uels {uelstr}];
% end
% gamsparams.n.name = 'n';
% gamsparams.n.uels = [];
% for jm = 1 : mop.od*2
%     uelstr = sprintf('%d',jm);
%     gamsparams.n.uels = [gamsparams.n.uels {uelstr}];
% end
% gamsparams.k.name = 'k';% number of objectives
% gamsparams.k.uels = [];
% for jm = 1 : mop.od
%     uelstr = sprintf('%d',jm);
%     gamsparams.k.uels = [gamsparams.k.uels {uelstr}];
% end
% 
% gamsparams.SumOneg.name = 'SumOne';
% gamsparams.SumOneg.val = mocpo_params.SumOne;
% gamsparams.SumOneg.type = 'parameter';
% 
% gamsparams.FloorC.name = 'FloorC';
% gamsparams.FloorC.val = least_multipliers;
% gamsparams.FloorC.type = 'parameter';
% 
% gamsparams.RoundLot.name = 'RoundLot';
% gamsparams.RoundLot.val = mocpo_params.Rl;
% gamsparams.RoundLot.type = 'parameter';

end


