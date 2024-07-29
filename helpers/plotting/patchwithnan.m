function patchwithnan(x,y1,y2,color,axhandle,varargin)
if nargin < 5 || isempty(axhandle)
    axhandle = gca; 
end
if nargin > 5
    extravars = reshape(varargin',2,[])'; 
    if ~any(cellfun(@(x) strcmp(x,'FaceAlpha'),extravars))
        extravars = [extravars; {'FaceAlpha',0.2}]; 
    end
    
    if ~any(cellfun(@(x) strcmp(x,'EdgeColor'),extravars))
        extravars = [extravars; {'EdgeColor','none'}]; 
    end
else
    extravars = []; 
end



y1 = y1(:)'; 
y2 = y2(:)'; 
x = x(:)'; 

datalocs = find(~isnan(y1)); 
dataends = [datalocs(diff(datalocs)>1), datalocs(find(datalocs,1,'last'))]; 
datastarts = [datalocs(find(datalocs,1,'first')), datalocs(find(diff(datalocs)>1)+1)]; 

regs = cell(1,length(datastarts)); 
for i = 1:(length(datastarts))
    regs{i} = datastarts(i):dataends(i); 
end


for i = 1:length(regs)
    
    xpatch = [x(regs{i}) fliplr(x(regs{i}))]; 
    ypatch = [y1(regs{i}), fliplr(y2(regs{i}))]; 

    if size(extravars,1)==1
        patch(axhandle,xpatch,ypatch,color,extravars{1,1},extravars{1,2});
    elseif size(extravars,1)==2
        patch(axhandle,xpatch,ypatch,color,extravars{1,1},extravars{1,2},extravars{2,1},extravars{2,2});
    elseif size(extravars,1)==3
        patch(axhandle,xpatch,ypatch,color,0.2,extravars{1,1},extravars{1,2},extravars{2,1},extravars{2,2},extravars{3,1},extravars{3,1});
    else
        patch(axhandle,xpatch,ypatch,color,'EdgeColor','none','FaceAlpha',0.2);
    end
end