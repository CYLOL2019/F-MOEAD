function init(mop)
%Set up the initial setting for the ICDV-MOEA/D.
%loading the params.
global params parDim objDim optimal_inds;  
  parDim = mop.pd;
  objDim = mop.od;
  %initial the subproblem's initital state.
  params.W = initweight(mop.od, params.popsize, strcmp(params.dmethod, 'te'));
  v        = squareform(pdist(params.W'));
  [~, params.neighbor]= sort(v);
  optimal_inds = randompoint(mop, params.popsize);
  optimal_inds = initidealpoint(optimal_inds);
end

%% initweight
function W = initweight(objDim, N, ifadjust)

U = floor(N^(1/(objDim-1)))-2;
M = 0;
while M<N
    U = U+1;
    M = noweight(U, 0, objDim); 
end

W   = zeros(objDim, M);
C   = 0;
V   = zeros(objDim,1);
[W, C] = setweight(W, C, V, U, 0, objDim, objDim);
W   = W / (U+0.0);

pos     = (W < 1.0E-5);
W(pos)  = 1.0E-5;

if ifadjust > 0
    W = 1.0./W;
    SW= sum(W);
    W = W./repmat(SW,objDim,1);
end

end

%%
function M = noweight(unit, sum, dim)

M = 0;

if dim == 1
    M = 1; 
    return;
end

for i=0:1:(unit - sum)
    M = M + noweight(unit, sum+i, dim-1);
end

end

%%
function [w, c] = setweight(w, c, v, unit, sum, objdim, dim)

if dim == objdim
    v = zeros(objdim, 1);
end

if dim == 1
    c       = c+1;
    v(1)    = unit-sum;
    w(:,c)  = v;
    return;
end

for i=0:1:(unit - sum)
    v(dim)  = i;
    [w, c]  = setweight(w, c, v, unit, sum+i, objdim, dim-1);
end

end