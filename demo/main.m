% Main function for F-MOEA/D
function main()
%clc
%clear
% SetPATH;
global mop params mocpo_params evalCounter optimal_inds end_marketime start_marketime T1 T2;
path('../public',path);
path('../public/NSGA2',path);
% path('../../../Finance_PyEnv',path);
problems = {'port1','port2','port3','port4','port5',...
    'port6','port7','port8','port9','port10',...
    'port11','port12','port13','port14','port15',...
    'port16','port17','port18','port19','port20'};
dirtime = '20090101,20201211_CVaR';
filename = sprintf('../%s/marketime.csv',dirtime);
marketime = load(filename);%marketime

nproblem = 10;
run_index =1;

for pn = 1 %1 : nproblem
    start_marketime = num2str(marketime(pn,1));
    end_marketime = num2str(marketime(pn,2));
    for r = run_index%1 : nrun
        % paramters for F-MOEAD
        params = loadparams_fmoead();    
        % parameters for the multiobjective constrained portfolio optimization
        [mocpo_params] = loadparams_mocpo(char(problems(1)),dirtime,pn);
        evalCounter = 0;
        tic;
        %init
        T1 = clock;
        init(mop)
        % main loop
        spareto = [optimal_inds.objectives];
%         ps = [optimal_inds.parameter;optimal_inds.combination];
%         pf = [optimal_inds.objectives];
%         archive = optimal_inds;
%         aps = ps;
%         apf = pf;
%         fes = evalCounter;
%         count_offsprings = params.popsize
        T2 = clock;
        while ~terminate()
%             evalCounter         
            [pareto, apareto] = fmoead(mop);
            spareto = [spareto;pareto.objectives];
            T2 = clock;
%             count_offsprings = count_offsprings + params.popsize
        end
        mocpo_params.grid = params.popsize;
        archive = evaluate_mocpo_final(pareto);
        pf = [archive.objectives];
        ps = [archive.parameter];
        
        sdir = sprintf("../data_2000s/port%d/run%d",pn,r);
        if ~exist(sdir, 'dir')
           mkdir(sdir)
        end
        runtime = toc;
        sfile = sprintf('%s/data.mat',sdir);
        disp(sfile);
        save(sfile,'pf','ps','spareto','runtime'); 
    end
       
end


end

function y =terminate()
    global params evalCounter count_offsprings T1 T2;
    y = etime(T2,T1)>= params.runtime;
%     y = evalCounter>=params.evaluation;
end

