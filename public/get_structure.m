function str = get_structure( name )
%STRUCTURE Summary of this function goes here
% 
% Structure used in this toolbox.
% 
% individual structure:
% parameter: the parameter space point of the individual. it's a column-wise
% vector.
% objective: the objective space point of the individual. it's column-wise
% vector. It only have value after evaluate function is called upon the
% individual.
%
% subproblem structure:
% weight: the decomposition weight for the subproblem.
% neighbour: the index of the neighbour of this subproblem.
% optimal: the current optimal value of the current structure. 
% curpoiont: the current individual of the subproblem.
% oldpoint: the backup individual of the subproblem used to caculate the
% utility.
% utility: the searching utility of this subproblem.
%
% testmop structure:
% name: the name of the test problem.
% od: the number of objective.
% pd: the number of variable.
% domain: the domain, which is a pd*2 matrix, with domain(:,1) and
% domain(:,2) to specificy the lower and upper limit on every variable.
%
% parameter strucutre.
% the parameter setting structure is explained in loadpapams.m
%
%

switch name
    case 'fmoead_parameter'
        str = struct('H',99,'dmethod', 'te','seed', 0,'Tm',10,'Tr',10,...
        'evaluation', 10^5, 'F',0.5, 'CR', 0.9,'neighbormating',0.9,'updatesize',2);
    case 'individual_MOCPO' 
        str = struct('parameter',[],'combination',[],'objectives',[],'estimation',[]);
    case 'individual_MOCPO_each' 
        str = struct('parameter',[],'combination',[],'weights',[],'objective',[],'estimation',[]);
    case 'individual_MultiScale' 
        str = struct('parameter',[],'combination',[],'objectives',[],'weights',[]);
    case 'individual' 
        str = struct('parameter',[],'objective',[],'estimation',[]);
    case 'subproblem' 
        str = struct('weight',[],'Tm', 10,'Tr', 10,'curpoint', [], ...
        'oldpoint',[],'optpoint',[],'mutsel',1,'fail', 0);
    case 'testmop'
        str = struct('name',[],'od',[],'pd',[],'domain',[],'func',[]);

    otherwise
        error('the structure name requried does not exist!');
end