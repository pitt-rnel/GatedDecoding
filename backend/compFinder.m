function Comp_info = compFinder(TS)

%%%% TRAINING SET
allfactors = cell2mat(cellfun(@(x) x.Factors,TS,'uni',0))'; % factors for training set
badgrasp = any(isnan(allfactors),2); Xtrain_grasp = allfactors(~badgrasp,:); % Find and remove NaNs

% Find click onsets and offsets
ptrace = cell2mat(cellfun(@(x) sum(x.Kin.Pos([7 10:15],:),1),TS,'uni',0));
ptrace = interp1(find(~isnan(ptrace)),ptrace(~isnan(ptrace)),1:length(ptrace));
ptrace = conv2(ptrace,ones(1,25)./25,'same');
pthresh = mean([prctile(unique(ptrace),95) prctile(unique(ptrace),5)]);
trlinds = cell2mat(cellfun(@(x) repmat(x.trial,1,size(x.Factors,2)),TS,'uni',0)); trlinds(badgrasp) = [];

offsettrace = ptrace(1:end-1)>pthresh & ptrace(2:end)<=pthresh;
onsettrace = ptrace(1:end-1)<pthresh & ptrace(2:end)>=pthresh;
if sum(onsettrace)~=sum(offsettrace)
    if find(onsettrace,1,'first')<find(offsettrace,1,'first') % start with click ON
        onsettrace(find(onsettrace,1,'last')) = false;
    else
        offsettrace(find(offsettrace,1,'first')) = false;
    end
end

clicktrace = zeros(1,length(ptrace)); clicktrace(onsettrace) = 1; clicktrace(offsettrace) = -1;
clicktrace(badgrasp) = [];

clickthreshs = [0.5 -0.5];
click_ons = find(clicktrace(1:end-1)<clickthreshs(1) & clicktrace(2:end)>clickthreshs(1) & diff(clicktrace)>0)';
click_offs = find(clicktrace(1:end-1)>clickthreshs(2) & clicktrace(2:end)<clickthreshs(2) & diff(clicktrace)<0)';

if length(click_ons)>length(click_offs)
    nmin = length(click_offs); 
    new_click_ons = NaN(nmin,1); 
    for i = 1:nmin
        new_click_ons(i) = click_ons(find(click_ons < click_offs(i),1,'last')); 
    end
    click_ons = new_click_ons; 
elseif length(click_ons)<length(click_offs)
    nmin = length(click_ons); 
    new_click_offs = NaN(nmin,1); 
    for i = 1:nmin
        new_click_offs(i) = click_ons(find(click_offs > click_ons(i),1,'first')); 
    end
    click_offs = new_click_offs; 
end
    
%%% Align components
graspdurs = cellfun(@(x) max([NaN,click_offs(find(click_offs>x,1,'first'))-x]),num2cell(click_ons));
restdurs = cellfun(@(x) max([NaN,click_ons(find(click_ons>x,1,'first'))-x]),num2cell(click_offs));
prepostpad = min(round(nanmedian(restdurs)./2),75);
durpad = min(round(nanmedian(graspdurs)/2),100);
onpad = -prepostpad:durpad;
offpad = -durpad:prepostpad;
facnum = size(Xtrain_grasp,2);
tlen = size(Xtrain_grasp,1);
onid = find(onpad==0);
offid = find(offpad==0)+length(onpad);
medid = round(nanmean([onid offid]));

click_on_reg = cellfun(@(x) onpad + x, num2cell(click_ons),'uni',0);
FacOnset = cellfun(@(x) [NaN(sum(x<1),facnum); Xtrain_grasp(x(x>0 & x<=tlen),:); NaN(sum(x>tlen),facnum)],...
    click_on_reg,'uni',0)';
FacOnsetAv = nanmean(cell2mat(reshape(FacOnset,1,1,[])),3);

click_off_reg = cellfun(@(x) offpad + x, num2cell(click_offs),'uni',0);
FacOffset = cellfun(@(x) [NaN(sum(x<1),size(Xtrain_grasp,2)); Xtrain_grasp(x(x>0 & x<=tlen),:); NaN(sum(x>tlen),facnum)],...
    click_off_reg,'uni',0)';

FacOffsetAv = nanmean(cell2mat(reshape(FacOffset,1,1,[])),3);

OnOffAv = [FacOnsetAv; FacOffsetAv];

SS = UniqueSubspaceSplitter({FacOnsetAv, FacOffsetAv},90,'do_plot',false); 

if size(SS.unique.C1,2)>1
    [R1] = varimax_sparsity(OnOffAv*SS.unique.C1,1:10); 
else
    R1 = 1; 
end
if size(SS.unique.C2,2)>1
    [R2] = varimax_sparsity(OnOffAv*SS.unique.C2,1:10); 
else
    R2 = 1; 
end

if size(SS.unique.C1,2)==0 || size(SS.unique.C2,2)==0
    [Rots,Comps] = varimax_sparsity(OnOffAv,1:10); 
else
    Rots = [SS.unique.C1*R1,SS.unique.C2*R2]; 
    Comps = OnOffAv*Rots; 
    Comps = Comps - repmat(nanmean(Comps(1:10,:)),size(Comps,1),1); 
end
    
SelfCross = deal(NaN(size(Comps,2),1));
TC = size(Comps,1);
T1 = 1:medid; T2 = (medid+1):TC;
for i = 1:size(Comps,2)
    c_i = Comps(:,i) - nanmax(Comps([1:10, (end-9):end],i));
    
    OnOffmag = [max(max(c_i(T1)),0)  max(max(c_i(T2)),0) ];    
    SelfCross(i) = -diff(OnOffmag);
end

on_comp = find(SelfCross==max(SelfCross));
off_comp = find(SelfCross==min(SelfCross));

trial_comps = cell2mat([reshape(cellfun(@(x) x*Rots(:,[on_comp,off_comp]),FacOnset,'uni',0),1,1,[]);...
    reshape(cellfun(@(x) x*Rots(:,[on_comp,off_comp]),FacOffset,'uni',0),1,1,[])]);

Comp_info.trial_comps = trial_comps;
Comp_info.Comps = Comps;
Comp_info.Rots = Rots;
Comp_info.on_comp = on_comp;
Comp_info.off_comp = off_comp;