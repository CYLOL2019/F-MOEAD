function [combination] = decode (population,PreAss,K)
[NoA, NP] = size(population);
[NULL, index] = sort(population,1,'descend');
combination = index(1:K,:);
for i = 1 : NP
    if ~ismember(PreAss,combination(:,i))
        combination(K,i) = PreAss;
    end
end
combination = sort(combination,1);
end
