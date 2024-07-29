function [metric,ldaplane,projs] = LDAmetric(X,L,numdisc)

if iscell(X)
    X = X(:); 
    Labcell = cell(length(X),1); 
    for i = 1:length(X)
        Labcell{i} = i*ones(length(X{i}),1); 
    end
    X = cell2mat(X); 
    L = cell2mat(Labcell); 
end

if nargin < 3
    numdisc = size(X,2)-1; 
end

g = ~any(isnan([X,L]),2); 
X = X(g,:); L = L(g); 

D = size(X,2); 
Lunique = unique(L); 
nL = length(Lunique); 


Cmeans = zeros(nL,D); 
for i = 1:nL
    Cmeans(i,:) = mean(X(L==Lunique(i),:)); 

end
totalmean = mean(Cmeans); 

% Between class covariances
Sb = cell(nL,1); 
for i = 1:nL
    Sb{i} = (Cmeans(i,:)-totalmean)'*(Cmeans(i,:)-totalmean); 
    Sb{i} = Sb{i} * sum(L==Lunique(i)); % multiply by # of class samples
end

SB = sum(cell2mat(reshape(Sb,1,1,[])),3); 

% Within class covariances
Sw = cell(nL,1); 
for i = 1:nL
    Ci = find(L==Lunique(i)); 
    swi = NaN(D,D,length(Ci)); 
    for j = 1:length(Ci)
        swi(:,:,j) = (X(Ci(j),:)-Cmeans(i,:))'*(X(Ci(j),:)-Cmeans(i,:)); 
    end
    Sw{i} = sum(swi,3); 
end

SW = sum(cell2mat(reshape(Sw,1,1,[])),3); 

%
[eigvecs,eigvals] = eig(SB,SW); 

eigvals = diag(eigvals); 

[~,ind] = sort(eigvals,'descend'); 

eigvecs = eigvecs(:,ind); 

metric = det(eigvecs(:,1:numdisc)' * SB * eigvecs(:,1:numdisc)) / ...
         det(eigvecs(:,1:numdisc)' * SW * eigvecs(:,1:numdisc)); 

ldaplane = eigvecs; 
projs = X*eigvecs; 
end



