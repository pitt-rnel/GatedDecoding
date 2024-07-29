SUBS = {'P2','P3'}; 
TransOccurring = Trans; 
for i = 1:length(Comps)
    TransOccurring{i} = cellfun(@(x) cellfun(@(y) y>0.05,x,'uni',0),Trans{i},'uni',0); 
end
nanmean50 = @(x) find(mean(isnan(x))<0.5); 
clrs = lines(10);  
[RVAR,ALIGNCOND] = deal(cell(length(Comps),1)); 
[compinfo,onsetalign,onsetsess] = deal(cell(length(Comps),4)); 
for z = 1:length(Comps)
    compnan = cellfun(@(x) {x{1}; NaN(20,size(x{1},2)); x{2}},Comps{z},'uni',0); 
    condvarnan = cellfun(@(x) {x{1}'; NaN(20,size(x{1},1)); x{2}'},Cvars{z},'uni',0); 
    transnan = cellfun(@(x) {x{1}'*1; NaN(20,size(x{1},1)); x{2}'*1},TransOccurring{z},'uni',0); 
    cat = cellfun(@(x) cell2mat(x),compnan,'uni',0); 
    cat_on = cellfun(@(x) x(:,1),cat,'uni',0); 
    cat_off = cellfun(@(x) x(:,2),cat,'uni',0); 
    
    mlen = max(cellfun(@length,cat_on)); 
    plenOn = max(cellfun(@(x) find(x==max(x),1,'first'),cat_on)); 
    plenOff = max(cellfun(@(x) length(x)-find(x==max(x),1,'first'),cat_off)); 

    [padded_on,padded_off,padded_on2,padded_off2] = deal(NaN(length(cat_on),3*mlen)); 
    [isOn,isOff] = deal(cell(length(cat_on),1)); 
    for i = 1:length(cat_on)
        pl = find(cat_on{i}==max(cat_on{i}),1,'first'); %peak location
        isOn{i} = (2*plenOn-pl+1):(2*plenOn-pl+length(cat_on{i})); 
        padded_on(i,isOn{i}) = cat_on{i}; 

        pl = find(cat_off{i}==max(cat_off{i}),1,'first'); %peak location
        isOff{i} = (3*mlen-pl+1):(3*mlen-pl + length(cat_off{i})); 
        padded_off(i,isOff{i}) = cat_off{i}; 
    end
    allnan = all(isnan(padded_on)); 
    padded_on(:,allnan) = []; 

    allnan = all(isnan(padded_off)); 
    padded_off(:,allnan) = []; 

    on_template = nanmean(padded_on); 
    off_template = nanmean(padded_off); 

    onpeak = find(on_template==max(on_template),1,'first'); 
    offpeak = find(off_template==max(off_template),1,'first'); 

    comp_columns = NaN(length(cat_on),2);
    [OnProjs,OffProjs,CondVars,TransHappen] = deal(cell(length(cat_on),2)); 
    for i = 1:length(cat_on)
        [OnProjs{i,1},OnProjs{i,2},OffProjs{i,1},OffProjs{i,2},CondVars{i,1},CondVars{i,2},TransHappen{i,1},TransHappen{i,2}] = deal(NaN(1,3*mlen)); 
        [rOn,rOff] = deal(zeros(size(cat{i},2),1)); 
        for j = 1:size(cat{i},2)
            g = ~any(isnan(compnan{i}{1}),2); 
            ti = (onpeak-floor(sum(g)/2)):(onpeak+floor(sum(g)/2)); ti = ti(1:sum(g)); 
            ti = ti - min([min(ti)-1, 0]); 
            rs = xcorr(compnan{i}{1}(g,j),on_template(ti),100,'normalized'); 
            rOn(j) = max(rs); 

            g = ~any(isnan(compnan{i}{3}),2); 
            ti = (offpeak-floor(sum(g)/2)):(offpeak+floor(sum(g)/2)); ti = ti(1:sum(g)); 
            ti = ti - min([min(ti)-1, 0]); 
            rs = xcorr(compnan{i}{3}(g,j),off_template(ti),100,'normalized'); 
            rOff(j) = max(rs); 
        end
        oncomp = find(rOn==max(rOn)); 
        offcomp = find(rOff==max(rOff)); 

        comp_columns(i,:) = [oncomp,offcomp]; 

        pl = find(compnan{i}{1}(:,oncomp)==max(compnan{i}{1}(:,oncomp)),1,'first'); %peak location
        is = (2*plenOn-pl+1):(2*plenOn-pl+size(compnan{i}{1},1)); 
        OnProjs{i,1}(1,is) = compnan{i}{1}(:,oncomp); 
        OffProjs{i,1}(1,is) = compnan{i}{1}(:,offcomp); 
        CondVars{i,1}(1,is) = condvarnan{i}{1}; 
        TransHappen{i,1}(1,is) = transnan{i}{1}; 

        pl = find(compnan{i}{3}(:,offcomp)==max(compnan{i}{3}(:,offcomp)),1,'first'); %peak location
        is = (2*mlen-pl+1):(2*mlen-pl + size(compnan{i}{3},1)); 
        OnProjs{i,2}(1,is) = compnan{i}{3}(:,oncomp); 
        OffProjs{i,2}(1,is) = compnan{i}{3}(:,offcomp); 
        CondVars{i,2}(1,is) = condvarnan{i}{3}; 
        TransHappen{i,2}(1,is) = transnan{i}{3}; 

    end
    nonanon = ~all(isnan( cell2mat(OnProjs(:,1)))); 
    nonanoff = ~all(isnan(cell2mat(OnProjs(:,2)))); 

    OnProjs(:,1) = cellfun(@(x) x(nonanon),OnProjs(:,1),'uni',0); 
    OnProjs(:,2) = cellfun(@(x) x(nonanoff),OnProjs(:,2),'uni',0); 

    OffProjs(:,1) = cellfun(@(x) x(nonanon),OffProjs(:,1),'uni',0); 
    OffProjs(:,2) = cellfun(@(x) x(nonanoff),OffProjs(:,2),'uni',0); 

    CondVars(:,1) = cellfun(@(x) x(nonanon),CondVars(:,1),'uni',0); 
    CondVars(:,2) = cellfun(@(x) x(nonanoff),CondVars(:,2),'uni',0); 

    TransHappen(:,1) = cellfun(@(x) x(nonanon),TransHappen(:,1),'uni',0); 
    TransHappen(:,2) = cellfun(@(x) x(nonanoff),TransHappen(:,2),'uni',0);      
    onav = nanmean(cell2mat(OnProjs(:,1))); onlocation = find(onav==max(onav),1,'first'); 
    offav = nanmean(cell2mat(OffProjs(:,2))); offlocation = find(offav==max(offav),1,'first'); 

    TransDuring = cellfun(@(x) any(x((onlocation+25):end)),TransHappen(:,1)) | cellfun(@(x) any(x(1:(offlocation-25))),TransHappen(:,2)); 

    % Force and digit variance
    force_indices = cellfun(@(x) ~isempty(x),ForConds{z}); 
    digit_indices = cellfun(@(x) ~isempty(x),DigitConds{z}); 

    f_trials = force_indices & ~TransDuring; 
    ft_trials = force_indices & TransDuring; 
    d_trials = digit_indices & ~TransDuring; 
    dt_trials = digit_indices & TransDuring; 

    trial_cell = {f_trials,ft_trials,d_trials,dt_trials}; 
    condnames = {'force','force+T','dig','dig+T'}; 

    x1 = nanmean50(cell2mat(OnProjs(:,1))); 
    x2 = nanmean50(cell2mat(OnProjs(:,2))); 
    x1p = x1; 
    x2p = x2 + x1(end) + 10; 

    [onpeak,infpeak,peakloc,offpeak,reg2_fromoff] = deal(NaN(1,4)); 
    [sp,xls,yls,RegVar] = deal(cell(2,4)); 
    figure('Position',[135 438 1579 502]); hold on; 
    for i = 1:length(trial_cell)
        if any(trial_cell{i})
            sp{1,i} = subplot(2,4,i); hold on; title([SUBS{z} ' ' condnames{i}]); 
            Y1on = bootstrap_LHM(@nanmean,cell2mat(OnProjs(trial_cell{i},1))); 
            Y2on = bootstrap_LHM(@nanmean,cell2mat(OnProjs(trial_cell{i},2))); 

            Y1off = bootstrap_LHM(@nanmean,cell2mat(OffProjs(trial_cell{i},1))); 
            Y2off = bootstrap_LHM(@nanmean,cell2mat(OffProjs(trial_cell{i},2))); 

            onpeak(i) = find(Y1on(:,3)==max(Y1on(:,3)),1,'first'); 
            offpeak(i) = find(Y2off(:,3)==max(Y2off(:,3)),1,'first'); 

            patchwithnan(x1p,Y1on(x1,1),Y1on(x1,2),clrs(1,:)); 
            patchwithnan(x2p,Y2on(x2,1),Y2on(x2,2),clrs(1,:)); 
            plot(x1p,Y1on(x1,3),'Color',clrs(1,:)); 
            plot(x2p,Y2on(x2,3),'Color',clrs(1,:));

            patchwithnan(x1p,Y1off(x1,1),Y1off(x1,2),clrs(2,:)); 
            patchwithnan(x2p,Y2off(x2,1),Y2off(x2,2),clrs(2,:)); 
            plot(x1p,Y1off(x1,3),'Color',clrs(2,:)); 
            plot(x2p,Y2off(x2,3),'Color',clrs(2,:));

            xls{1,i} = xlim; yls{1,i} = ylim; 

            sp{2,i} = subplot(2,4,i+4); hold on; 

            Y1 = bootstrap_LHM(@nanmean,cell2mat(CondVars(trial_cell{i},1))); 
            Y2 = bootstrap_LHM(@nanmean,cell2mat(CondVars(trial_cell{i},2))); 

            infpeak(i) = find(Y1(:,3)==max(Y1(:,3)),1,'first'); 

            patchwithnan(x1p,Y1(x1,1),Y1(x1,2),clrs(4,:)); 
            patchwithnan(x2p,Y2(x2,1),Y2(x2,2),clrs(4,:)); 
            plot(x1p,Y1(x1,3),'Color',clrs(4,:)); 
            plot(x2p,Y2(x2,3),'Color',clrs(4,:)); 

            xls{2,i} = xlim; yls{2,i} = ylim; 

            reg1 = onpeak(i) + 10 + (-2:2); 
 
            condvars1 = cell2mat(CondVars(trial_cell{i},1)); 
            condvars2 = cell2mat(CondVars(trial_cell{i},2)); 
            RegVar{1,i} = nanmean(condvars1(:,reg1),2); 
            RegVar{1,i}(RegVar{1,i}==0) = NaN; 

            reg2 = offpeak(i)-10+ (-2:2); 

            RegVar{2,i} = nanmean(condvars2(:,reg2),2); 

            reg3 = onpeak(i) + (-50:50); 

            plot(mean(reg1)*[1 1],yls{2,i},'r'); 
            plot(mean(reg2)*[1 1] + x1(end)+10,yls{2,i},'r'); 

            cv1_norm = (condvars1(:,reg3)-(repmat(min(condvars1(:,reg3),[],2),1,length(reg3))))./repmat(range(condvars1(:,reg3),2),1,length(reg3)); 
            compinfo{z,i} = cv1_norm; 
            onsetalign{z,i} = (Y1on(reg3,3)-min(Y1on(reg3,3)))./range(Y1on(reg3,3)); 

            onsetall = cell2mat(OnProjs(trial_cell{i},1)); 
            onsetsess{z,i} = onsetall(:,reg3); 
        else
            compinfo{z,i} = NaN(1,101); 
            onsetalign{z,i} = NaN(101,1); 
            onsetsess{z,i} = NaN(1,101); 
        end
    end

    xl = [min(cell2mat(xls(1,:))), max(cell2mat(xls(1,:)))]; 
    yl_1 = [min(cell2mat(yls(1,:))), max(cell2mat(yls(1,:)))];    
    yl_2 = [min(cell2mat(yls(2,:))), max(cell2mat(yls(2,:)))];

    for i = 1:length(trial_cell)
        if ~isempty(sp{1,i})
            set(sp{1,i},'XLim',xl,'YLim',yl_1); 
            set(sp{2,i},'XLim',xl); 
        end
    end

    make_pretty; 

    RVAR{z} = RegVar; 
end