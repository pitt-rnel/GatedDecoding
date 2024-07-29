%%
figure('Position',[53 149 1816 826]); hold on; 
mrks = {'.','o'}; 
info_reduction = cell(1,length(Comps)); 
P_vals = NaN(1,4); 
for i = 1:size(RVAR{1},2)
    subplot(1,size(RVAR{1},2),i); hold on; title(condnames{i}); 
    for z = [1 2]
        if ~isempty(RVAR{z}{1,i})
            plot(RVAR{z}{1,i},RVAR{z}{2,i},mrks{z},'Color','k'); 
            info_reduction{z}{i}([1 3 2]) = bootstrap_LHM(@(x) (x(~any(isnan(x),2),1)\x(~any(isnan(x),2),2))*100,[RVAR{z}{1,i},RVAR{z}{2,i}]); 
        else
            info_reduction{z}{i} = []; 
        end
    end
    xlabel('reg 1 cond var'); 
    ylabel('reg 2 cond var'); 
    xl = xlim; yl = ylim; zl = [0,max([xl yl])]; 

    xlim(zl); ylim(zl); 

    plot(zl,zl,'k:'); axis square;  

    if i==1; legend({'P2','P3'}); end

    [~,P_vals(1,i)] = ttest2(cell2mat(cellfun(@(x) x{1,i},RVAR([1 2]),'uni',0)), cell2mat(cellfun(@(x) x{2,i},RVAR([1 2]),'uni',0)));
end
make_pretty; 