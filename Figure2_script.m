loadPath = 'C:\Users\dekleva\OneDrive - University of Pittsburgh\Desktop\DataCondensed\GatedDecoding\'; 
subjectId = {'P2','P3'}; 
dataset_nums = {1:107, 1:63}; 

[Comps,Perfs,Cvars,Metrics,Trans,ForConds,DigitConds,TransConds] = deal(cell(2,1)); 
for i = 1:length(subjectId)

    [Comps{i},Perfs{i},Cvars{i},Metrics{i},Trans{i},ForConds{i},DigitConds{i},TransConds{i}] = ...
        ExtractGraspDatasetsComponents(loadPath,subjectId{i},dataset_nums{i});

end

%%
PlotTransientsAndVariance;

%%
PlotScatterVariance;

%%
PlotLags;