function testmop( testname, pdim)
%Get test multi-objective problems from a given name. 
%   The method get testing or benchmark problems for Multi-Objective
%   Optimization. The test problem will be encapsulated in a structure,
%   which can be obtained by function get_structure('testmop'). 
%   User get the corresponding test problem, which is an instance of class
%   mop, by passing the problem name and optional dimension parameters.
global mop;
    switch testname
        case {'ZDT1','ZDT2','ZDT3','ZDT4','ZDT5','ZDT6','test'}
            mop=ZDT(testname,pdim);
        case {'DTLZ1','DTLZ2','DTLZ3','DTLZ4','DTLZ5','DTLZ6','DTLZ7'}
            odim=3;
            mop=DTLZ(testname,pdim,odim);
        case {'zzj_f1','zzj_f2','zzj_f3','zzj_f4','zzj_f5','zzj_f6','zzj_f7','zzj_f8','zzj_f9','zzj_f10'}
            mop=zzj( testname, pdim);        
        case {'tec09_f1','tec09_f2','tec09_f3','tec09_f4','tec09_f5','tec09_f6','tec09_f7','tec09_f8','tec09_f9'}
            mop =tec09( testname, pdim );
        case {'uf1', 'uf2','uf3','uf4','uf5','uf6','uf7'}
            mop=cecproblems(mop, testname, pdim);
            mop.od=2;    
        case {'uf8','uf9','uf10'}
            mop=cecproblems(mop, testname, pdim);
            mop.od=3;   
        case {'r2_dtlz2_m5', 'r3_dtlz3_m5', 'wfg1_m5'}
            mop=cecproblems2(mop, testname, pdim); 
        case {'MOP1','MOP2','MOP3','MOP4','MOP5','MOP6','MOP7'}
            MOP(testname,pdim);
        case {'port1','port2','port3','port4','port5','port6','port7','port8','port9','port10',...
    'port11','port12','port13','port14','port15','port16','port17','port18','port19','port20'}
            mop = portfolio(testname, pdim);
        otherwise 
            error('Undefined test problem name'); 
    end    
end
% Portfolio Optimization D1-D20
function p= portfolio(testname,dim)
 p.name=upper(testname);
 p.pd=dim;
 p.od=2;
 p.domain = [zeros(dim,1),ones(dim,1)];

 p.func = @portpfeval;
    % portfolio evaluation function
    % x and y are columnwise, the imput x must be inside the search space and
    % it could be a matrix
    function y = portpfeval(x)
        global u Covariance
        [dim, num]  = size(x);
        y = zeros(2,num);
        y(1,:) = sum(sum(x*x'.*Covariance));
        y(2,:) = -x'*u;
%         for i = 1 : num
%             y(1,i) = bsxfun(funcov, x(:,i), Covariance);
%         end

    end
end
%cec09 UF1 - UF10
function p=cecproblems(p, testname,dim)
 p.name=upper(testname);
 p.pd=dim;
 
 p.domain=xboundary(upper(testname),dim);
 %p.domain = [zeros(dim,1),ones(dim,1)];
 p.func=cec09(upper(testname));
end

%cec09 UF11 - UF13
function p=cecproblems2(p, testname,dim)
 p.name=upper(testname);
 p.pd=dim;
 p.od=2;
 
 p.domain=xboundary(upper(testname),dim);
 %p.domain = [zeros(dim,1),ones(dim,1)];
 p.func=cec09m(upper(testname));
end

