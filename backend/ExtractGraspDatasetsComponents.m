function[Comps,PerfSplit,CvarSplit,MetricSplit,TransSplit,ForceConditions,DigitConditions,TransConditions] = ExtractGraspDatasetsComponents(data_dir,subID,dataset_range)

N = length(dataset_range); 
[Comps,PerfSplit,CvarSplit,MetricSplit,TransSplit,ForceConditions,DigitConditions,TransConditions] = deal(cell(N,1)); 
data_dir = [data_dir subID filesep]; 
listing = dir(data_dir); 
for i = 1:N

    clc; fprintf('%d/%d\n',i,N); 
    % try
        fili = dataset_range(i); 
        matching_file_loc = cellfun(@(x) contains(x,[subID '_Data_' num2str(fili) '_']),{listing.name}); 
        load([data_dir listing(matching_file_loc).name],'TS'); 
              
        Comp_info = compFinder(TS); % Find transients
        
        forcesinc = unique(cell2mat(cellfun(@(x) x.Force_levels(x.Force_levels>0),TS,'uni',0))); 
        targs = cell2mat(cellfun(@(x) x.target,TS','uni',0)); 
        fingersinc = find(sum(targs(:,10:14),'omitnan')>0); 
        transinc = find(sum(targs(:,1:3),'omitnan')>0); 
        ForceConditions{i} = forcesinc; 
        DigitConditions{i} = fingersinc; 
        TransConditions{i} = transinc; 

        OnOff = findOnsetOffset(TS); 
        [FacSplit,FacCat,Ksplit] = splitFacsOnOff(TS,OnOff); 
        binsizes = cellfun(@(x) size(x,2),FacSplit(1,:)); 
        SpdSplit = cellfun(@(x) sqrt(sum(x.^2,1)),Ksplit,'uni',0); 
        spd = SpeedByTime(SpdSplit);

        avfac = mean(cell2mat(reshape(FacCat,1,1,[])),3,'omitnan')'; 
        avfac = avfac-mean(avfac(1:5,:),'omitnan'); 

        Ssparse = avfac*Comp_info.Rots; 
        OnOff = [Comp_info.on_comp,Comp_info.off_comp]; 
        OtherComps = 1:size(Ssparse,2); OtherComps(ismember(OtherComps,OnOff)) = []; 
        Ssparse = Ssparse(:,[OnOff OtherComps]); 
        
        Comps{i} = mat2cell(Ssparse,[sum(binsizes(1:2)),sum(binsizes(3:4))],size(Ssparse,2)); 

        if ~isempty(forcesinc)

            [perf,cvar,metric] = ForceByTime(TS,FacSplit);

        end

        if ~isempty(fingersinc)
            for ii = 1:length(TS)
                finger_disc = nansum(TS{ii}.target(:,10:14))>0; 
                fingers = find(finger_disc);
                
                TS{ii}.digits = fingers; 
                TS{ii}.digit_id = str2double(sprintf('%d',fingers)); 
            end
            [perf,cvar,metric] = DigitByTime(TS,FacSplit);
        end

        PerfSplit{i} = {cell2mat(perf(1:2)), cell2mat(perf(3:4))}; 
        CvarSplit{i} = {cell2mat(cvar(1:2)), cell2mat(cvar(3:4))}; 
        MetricSplit{i} = {cell2mat(metric(1:2)), cell2mat(metric(3:4))}; 
        TransSplit{i} = {cell2mat(spd(1:2)), cell2mat(spd(3:4))}; 
    % catch
    %     fprintf('something wrong on set %d\n',i); 
    % end

end

end

function[OnOff] = findOnsetOffset(TS)

    getOn = @(x,th) min([find((x(1:(end-1)) < th) & (x(2:end) >= th))+1, NaN]); 
    getOff = @(x,th) min([find((x(1:(end-1)) > th) & (x(2:end) <= th)),NaN]); 

    prange = [min(cell2mat(cellfun(@(x) reshape(sum(x.Kin.Pos([7 10:14],:)),1,[]),TS,'uni',0))), ...
              max(cell2mat(cellfun(@(x) reshape(sum(x.Kin.Pos([7 10:14],:)),1,[]),TS,'uni',0))) ]; 

    thresh = median(prange); 

    OnOff = cell2mat(cellfun(@(x) [getOn(sum(x.Kin.Pos([7 10:14],:)),thresh),getOff(sum(x.Kin.Pos([7 10:14],:)),thresh)],TS,'uni',0)'); 
end

function[FacSplit,FacCat,Ksplit,Kcat] = splitFacsOnOff(TS,OnOff)
        n = size(OnOff,1); 
        hold_lens = OnOff(:,2)-OnOff(:,1); 
        minhold = min(hold_lens); 
        minpad1 = min(OnOff(:,1))-1; 
        minpad2 = min(cellfun(@(x) size(x.Factors,2),TS)' - OnOff(:,2))-1; 

        prelen = min(minpad1,150); 
        holdlen = floor(minhold/2); 
        postlen = min(minpad2,150); 

        regs = {repmat(-prelen:-1,n,1) + OnOff(:,1),...
                repmat(0:holdlen,n,1) + OnOff(:,1),...
                repmat(-holdlen:-1,n,1) + OnOff(:,2),...
                repmat(0:postlen,n,1) + OnOff(:,2)}; 

        [FacSplit,Ksplit] = deal(cell(length(TS),length(regs))); 
        [FacCat,Kcat] = deal(cell(length(TS),1)); 
        for i = 1:length(TS)
            for j = 1:length(regs)
                FacSplit{i,j} = TS{i}.Factors(:,regs{j}(i,:)); 
                Ksplit{i,j} = TS{i}.Kin.Vel(1:3,regs{j}(i,:)); 
            end
            FacCat{i} = cell2mat(FacSplit(i,:)); 
            Kcat{i} = cell2mat(Ksplit(i,:)); 
        end

end

function[perf,cvar,metric] = ForceByTime(TS,FacSplit)

    L = cellfun(@(x) x.Force_levels(x.Force_levels>0),TS)'; 
    uL = unique(L); 
    nd = length(uL)-1; 

    [metric,perf,cvar] = deal(cell(1,size(FacSplit,2))); 
    for i = 1:size(FacSplit,2) %epoch
        for j = 1:size(FacSplit{1,i},2) %timepoint 
            X = cell2mat(cellfun(@(x) x(:,j)',FacSplit(:,i),'uni',0)); 
            [metric{i}(j),~,~] = LDAmetric(X,L,nd);
            [~,~,perf{i}(j)] = LDA_LOO(X,L); 

            Xc = cell2mat(cellfun(@(x) mean(X(L==x,:)),num2cell(uL),'uni',0)); 
            cvar{i}(j) = sum(var(Xc,[],1)); 

        end
    end
end

function[perf,cvar,metric] = DigitByTime(TS,FacSplit)

    L = cellfun(@(x) x.digit_id,TS)'; 
    uL = unique(L); 
    nd = length(uL)-1; 

    [metric,perf,cvar] = deal(cell(1,size(FacSplit,2))); 
    for i = 1:size(FacSplit,2) %epoch
        for j = 1:size(FacSplit{1,i},2) %timepoint 
            X = cell2mat(cellfun(@(x) x(:,j)',FacSplit(:,i),'uni',0)); 
            [metric{i}(j),~,~] = LDAmetric(X,L,nd);
            [~,~,perf{i}(j)] = LDA_LOO(X,L); 

            Xc = cell2mat(cellfun(@(x) mean(X(L==x,:)),num2cell(uL),'uni',0)); 
            cvar{i}(j) = sum(var(Xc,[],1)); 

        end
    end
end

function[spd] = SpeedByTime(SpdSplit)

    [spd] = deal(cell(1,size(SpdSplit,2))); 
    for i = 1:size(SpdSplit,2) %epoch
        for j = 1:size(SpdSplit{1,i},2) %timepoint 
            spd{i}(j) = mean(cellfun(@(x) x(:,j)',SpdSplit(:,i))); 
        end
    end
end


function[dateOut] = extractDate(Data)
    
    fstr = Data.Files{1}; 
    datedashes = strfind(fstr,'-'); 
    dotlocs = strfind(fstr,'.'); 
    undlocs = strfind(fstr,'_'); 

    dstart = dotlocs(find(dotlocs<datedashes(1),1,'last'))+1; 
    dend = undlocs(find(undlocs>datedashes(2),1,'first'))-1; 

    dateOut = fstr(dstart:dend); 
end
