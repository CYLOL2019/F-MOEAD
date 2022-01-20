function obj = subobjective(weight, ind, idealpoint, method)
%SUBOBJECTIVE function evaluate a point's objective with a given method of
%decomposition. 

%   Two method are implemented by far is Weighted-Sum and Tchebesheff.
%   weight: is the decomposition weight.(column wise vector).
%   ind: is the individual point(column wise vector).
%   idealpoint: the idealpoint for Tchebesheff decomposition.
%   method: is the decomposition method, the default is 'te' when is
%   omitted.
%   
%   weight and ind can also be matrix. in which have two scenairos:
%   When weight is a matrix, then it's treated as a column wise set of
%   weights. in that case, if ind is a size 1 column vector, then the
%   subobjective is computed with every weight and the ind; if ind is also
%   a matrix of the same size as weight, then the subobjective is computed
%   in a column-to-column, with each column of weight computed against the
%   corresponding column of ind. 
%   A row vector of subobjective is return in both case.

% idealpoint=idealpoint*(1-0.5^(itrCounter*0.5));
    if (nargin==2)
        obj = ws(weight, ind);
    elseif (nargin==3)
        obj = te(weight, ind, idealpoint);
    else
        if strcmp(method, 'ws')
            obj=ws(weight, ind);
        elseif strcmp(method, 'te')
            obj=te(weight, ind, idealpoint);
        elseif strcmp(method, 'pbi')
            obj=pbi(weight, ind, idealpoint);
        elseif strcmp(method,'ww')
            obj=ww(weight,ind);
        else
            obj= te(weight, ind, idealpoint);
        end
    end
end

% weighted sum scalarization Function
function obj = ws(weight, ind)
    obj =sum(weight.*ind,1);
end

% Techbycheff Scalarization Function%%%   ind--目标值  weight--权值
function obj = te(weight, ind, idealpoint)
global mocpo_params nadirpoint
    s = size(weight, 2);     %权值个数
    indsize = size(ind,2);   %目标个数
    weight((weight == 0))=0.00001;
    ind = (ind - idealpoint)./(nadirpoint - idealpoint);
    %calculation for subproblems when objectives are PF
    if indsize == mocpo_params.grid
        part2 = abs(ind);
        obj = zeros(1,s);
        for i = 1 : s
            %First, max means te for each points on the obtained PF
            %Second, min means the minimum of te for one subset
            obj(i) = min(max(weight(:,i).*part2));
        end
    elseif indsize == s* mocpo_params.grid
        obj = zeros(1,s);
        for i = 1 : s
            % (i-1)*subpopsize:i*subpopsize is the index of i-th PF
            part2 = abs(ind(:, (i-1)*mocpo_params.grid+1:i*mocpo_params.grid));                
            %First, max means te for each points on the obtained PF
            %Second, min means the minimum of te for one subset
            obj(i) = min(max(weight(:,i).*part2));
        end
    end

    
%     if indsize==s 
%         part2 = abs(ind-idealpoint(:,ones(1, indsize)));
%         obj = max(weight.*part2);
%     elseif indsize ==1
%         part2 = abs(ind-idealpoint);
%         obj = max(weight.*part2(:,ones(1, s)));
%     elseif s==1
%         part2=abs(ind-idealpoint(:,ones(1, indsize)));
%         obj = max(weight(:,ones(1, indsize)).*part2);
%         
%     else
%         error('individual size must be same as weight size, or equals 1');
%     end
    
end
function obj=pbi(weight,ind,idealpoint)

    s=size(weight,2);
    indsize=size(ind,2);
  weight((weight==0))=0.0000001;
    sita=3;

    if indsize==s
        d1=abs(sum((ind-idealpoint(:,ones(1,indsize))).*weight,1))./(sqrt(sum(weight.*weight,1))); 
        dd=ind-(idealpoint(:,ones(1,s))+d1(ones(size(weight,1),1),:).*weight);
        d2=sqrt(sum(dd.*dd,1));
        obj=d1+sita*d2;
    elseif indsize==1
        part2=ind-idealpoint;
        d1=abs(sum(part2(:,ones(1,s)).*weight,1))./sqrt(sum(weight.*weight,1));
        dd=ind(:,ones(1,s))-(idealpoint(:,ones(1,s))+d1(ones(size(weight,1),1),:).*weight);
        d2=sqrt(sum(dd.*dd,1));
        obj=d1+sita*d2;
    elseif s==1
        d1=abs(sum((ind-idealpoint(:,ones(1,indsize))).*weight(:,ones(1, indsize)),1))./(sqrt(sum(weight(:,ones(1, indsize)).*weight(:,ones(1, indsize)),1))); 
        dd=ind-(idealpoint(:,ones(1,indsize))+d1(ones(size(weight,1),1),:).*weight(:,ones(1, indsize)));
        d2=sqrt(sum(dd.*dd,1));
        obj=d1+sita*d2;
    else
        error('individual size must be same as weight size, or equals 1');
    end   

end
function obj = ww(weight, ind)
obj=sum(weight.*(ind.^2),1);
end