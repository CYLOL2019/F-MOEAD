function ind = genetic_op_DEPM_swap(index,domain,useneighbour)%domin 是最新父代种群大小的意思
%GENETICOP function implemented the DE operation to generate a new
%individual from a subproblems and its neighbours.

%   subproblems: is all the subproblems.
%   index: the index of the subproblem need to handle.
%   domain: the domain of the origional multiobjective problem.
%   ind: is an individual structure.
  global params parDim mocpo_params optimal_inds;
%% %%% %%%%%%%%%
operator_sign = randi(2);
ind = get_structure('individual_MOCPO');
if operator_sign == 1
    % mutsel=subproblems(index).mutsel;
    mutsel=1;
    %%%%%%%参数假设%%%%%
      parents = mateselection(index,3,useneighbour);
        switch mutsel
            case 1
              parents = parents(1:3);
            case 3
               parents = parents(1:4);
            otherwise
       end
      selp=ceil(length(params.F)*rand);
      newpoint = de_crossover(parents, params.F(selp), params.CR(selp), domain,mutsel);

      ind.parameter = newpoint;
          % 
       while all(ind.parameter==optimal_inds(index).parameter)
        ind = realmutate(ind, domain, 1/parDim);
       end
        ind = realmutate(ind, domain, 1/parDim);
    %    ind = gaussian_mutate(ind, 1/parDim, domain);
% 
else

     %Swap
    swap_sign = randi(3);
    newpoint = optimal_inds(index).parameter;
    pretemp = newpoint;%set pre-assignment to 1
    pretemp(:,mocpo_params.PreAss) = 1;
    [NULL, index] = sort(pretemp,'descend'); 
    pindex = index(1:mocpo_params.K);
    while 1
        first = pindex(randi(mocpo_params.K));
        if isempty(intersect(first, mocpo_params.PreAss))
          break;
        end
    end
    %     [NULL first] = max(ind(1:end-2));

    switch (swap_sign)
        case 1
            %Highest return

           for traverseindex = 1 : parDim
               if ~ismember(mocpo_params.u_idx(traverseindex),pindex)
                   I = mocpo_params.u_idx(traverseindex);
                   break;
               end
           end
       case 2
            %Lowest Epexceted loss of t time invervals (related to CVaR)

           for traverseindex = 1 : parDim
               if ~ismember(mocpo_params.cvar_idx(traverseindex),pindex)
                   I = mocpo_params.cvar_idx(traverseindex);
                   break;
               end
           end
        case 3
           %Lowest loss at time T (related to VaR)

           for traverseindex = 1 : parDim
               if ~ismember(mocpo_params.var_idx(traverseindex),pindex)
                   I = mocpo_params.var_idx(traverseindex);
                   break;
               end
           end 


    end
    temp     = newpoint(first);
    newpoint(first) = newpoint(I);
    newpoint(I) = temp;
    ind.parameter = newpoint;
end
  


    % decode selection for GAMS
    ind.combination = decode(ind.parameter,mocpo_params.PreAss,mocpo_params.K);
  clear points selectpoints oldpoint randomarray deselect newpoint parentindex si;
end

function ind = realmutate(ind, domains, rate)
%REALMUTATE Summary of this function goes here
%   Detailed explanation goes here

  % double rnd, delta1, delta2, mut_pow, deltaq;
  % double y, yl, yu, val, xy;
  % double eta_m = id_mu;
  global rnduni parDim;

  eta_m=20;
  
  if (isstruct(ind))
      a = ind.parameter;
  else
      a = ind;
  end
  
  %[r rnduni] = crandom(rnduni);
  %id_rnd = ceil(r*parDim);
  
  for j = 1:parDim
      %[r rnduni] = crandom(rnduni);
      r = rand;
      
      if (r <= rate) 
        y = a(j);
        yl = domains(j,1);
        yu = domains(j,2);
        delta1 = (y - yl) / (yu - yl);
        delta2 = (yu - y) / (yu - yl);

        %[rnd rnduni] = crandom(rnduni);
        rnd = rand;
        mut_pow = 1.0 / (eta_m + 1.0);
        if (rnd <= 0.5) 
	      xy = 1.0 - delta1;
	      val = 2.0 * rnd + (1.0 - 2.0 * rnd) * (xy^(eta_m + 1.0));
	      deltaq = (val^mut_pow) - 1.0;
        else 
	      xy = 1.0 - delta2;
	      val = 2.0 * (1.0 - rnd) + 2.0 * (rnd - 0.5) * (xy^ (eta_m + 1.0));
	      deltaq = 1.0 - (val^mut_pow);
        end
	  
        y = y + deltaq * (yu - yl);
        if (y < yl)
	      y = yl;
        end
        if (y > yu)
	      y = yu;
        end
        
        a(j) = y;        
      end
  end
  
  if isstruct(ind)
      ind.parameter = a;
  else
      ind = a;
  end
end

function ind = gaussian_mutate( ind, prob, domain)
%GAUSSIAN_MUTATE Summary of this function goes here
%   Detailed explanation goes here

  if isstruct(ind)
      x = ind.parameter;
  else
      x  = ind;
  end

    parDim = length(x);
    lowend  = domain(:,1);
    highend =domain(:,2);
    sigma = (highend-lowend)./20;
    
    newparam = min(max(normrnd(x, sigma), lowend), highend);
    C = rand(parDim, 1)<prob;
    x(C) = newparam(C);
    
  if isstruct(ind)
      ind.parameter = x;
  else
      ind = x;
  end
end

% select the candidate for evoluationary mating, i.e. finding the DE parent.
function select = mateselection(index,size,useneighbour)

global params;

if useneighbour
  neighborhood = params.neighbor(1:params.Tm,index);
else
  neighborhood = params.neighbor(:,index);
end
select     = ones(size,1)*index;
while select(2)==select(1) || select(3)==select(1) || select(3)==select(2)
    select(2:3) = randsample(neighborhood, 2);
end
   
end

function ind = de_crossover(parents, F, CR, domain,mutsel)
  global  parDim optimal_inds;
  
  %  [r rnduni]=crandom(rnduni);
  r = rand;
  jrandom = ceil(r*parDim);
  
  %retrieve the individuals.
  points = [optimal_inds(parents)];
  selectpoints = [points.parameter];

  switch mutsel
        case 1
          cross = selectpoints(:,1) + F.*(selectpoints(:,3)-selectpoints(:,2));
        case 2
          cross = selectpoints(:,1) + F.*(selectpoints(:,3)-selectpoints(:,2))+F.*(selectpoints(:,5)-selectpoints(:,4));
        case 3
          cross = selectpoints(:,1) + F.*(selectpoints(:,3)-selectpoints(:,2))+rand.*(selectpoints(:,4)-selectpoints(:,1));
% cross = selectpoints(:,1) + 0.5*F.*(selectpoints(:,1)-selectpoints(:,3)-selectpoints(:,2));
        otherwise
   end

  ind = selectpoints(:,1);
  
  %DE operation.
  for i=1:parDim
      %[r rnduni]=crandom(rnduni);
      r = rand;
      if (r<CR || i==jrandom)
          ind(i) = cross(i);
      end
      
      %handle the boundary.
      lowbound = domain(i,1);
      upbound = domain(i,2);
      if (ind(i)<lowbound)
          %[r rnduni]=crandom(rnduni);
%           r = rand;
          ind(i)=lowbound+ r*(selectpoints(i,1)-lowbound);
%  ind(i)=lowbound;
      elseif (ind(i)>upbound)
          %[r rnduni]=crandom(rnduni);
%           r = rand;
          ind(i)=upbound - r*(upbound - selectpoints(i,1));
%           ind(i)=upbound ;
      end
  end
end