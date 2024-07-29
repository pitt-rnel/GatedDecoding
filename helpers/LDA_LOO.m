function [classout,classout_cell,perf,probout] = LDA_LOO(Xcell,Lcell)

if ~iscell(Xcell)
    Xcell = mat2cell(Xcell,ones(size(Xcell,1),1),size(Xcell,2)); 
    Lcell = num2cell(Lcell(:)); 
end
unl = unique(cell2mat(Lcell)); 
[classout_cell,ps_cell] = deal(cell(length(Xcell),1)); 
for i = 1:length(Xcell) % trials
    
    traini = 1:length(Xcell); traini(i) = []; 
    testi = i; 
    
    X = cell2mat(Xcell(traini)); 
    L = cell2mat(Lcell(traini)); 
    Xt = Xcell{testi}; 
    
    good = ~any(isnan([X,L]),2); 
    goodt = ~any(isnan(Xt),2); 
    
    [~,ps_c,C] = LDA_classify(X(good,:),L(good),[],Xt(goodt,:)); 
    unL = unique(L); for j = 1:length(unL); C(C==j)=unL(j); end
    
    classout_cell{i} = NaN(size(Xt,1),1); 
    classout_cell{i}(goodt) = C; 
%     ps_cell{i}(goodt,:) = zeros(1,length(goodt)); 
    ps_cell{i}(goodt,ismember(unl,L(good))) = ps_c; 
end
    
classout = cell2mat(classout_cell); 

perf = nanmean(classout==cell2mat(Lcell)); 
probout = cell2mat(ps_cell); 