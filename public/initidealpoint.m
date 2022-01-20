function [x] = initidealpoint(x)
global cachepoints idealpoint nadirpoint
% v = [];
iv = ones(2,1)*inf;
nv = ones(2,1)*-inf;
% get ideal risk and nadir return
for i  = 1 : length(x)
    vi   = cvar_weight(x(i).combination,[1;0]);
    if vi(1) < iv(1)
        iv(1) = vi(1);
        nv(2) = vi(2);
    end
end
% get ideal return and nadir risk
for i  = 1 : length(x)
    vi   = cvar_weight(x(i).combination,[0;1]);
    if vi(2) < iv(2)
        iv(2) = vi(2);
        nv(1) = vi(1);
    end
end
idealpoint = iv;
nadirpoint = nv;
cachepoints = [iv(1),nv(1);nv(2),iv(2)];%[ideal_risk   nadir_risk
                                        % nadir_return ideal_return]
for i  = 1 : length(x)
    vi   = cvar(x(i).combination);
    x(i).objectives = vi;

end
end
