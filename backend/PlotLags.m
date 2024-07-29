%%
clrs = lines(10); 
peaklag = cell(size(onsetsess)); 
for i = 1:numel(onsetsess)
    for j = 1:size(onsetsess{i},1)
        x1 = compinfo{i}(j,:)';
        x2 = onsetsess{i}(j,:)'; 
        g = ~any(isnan([x1 x2]),2); 
        [r,lag] = xcorr(x1(g),x2(g),50);
        peaklag{i}(1,j) = lag(find(r==max(r),1,'first')); 
    end
end
%%
pl1 = cell2mat(peaklag(1,:))'; 
pl2 = cell2mat(peaklag(2,:))'; 

figure; hold on; 
[x1,n1] = linehist(-50:5:50,pl1,'suppress_plot'); 
[x2,n2] = linehist(-50:5:50,pl2,'suppress_plot'); 

plot(x1,n1,'Color',clrs(4,:),'LineWidth',2); 
plot(x2,n2,':','Color',clrs(4,:),'LineWidth',2); 
yl = ylim; 
plot(median(pl1),yl(2),'k.'); 
plot(median(pl2),yl(2),'ko'); 
plot([0 0],yl,'k:'); 
xlim([-50 50]); 
xlabel('optimal onset-variance lag')
ylabel('# sessions'); 
title(sprintf('P2: %dms  P3: %dms',median(pl1)*20,median(pl2)*20)); 
make_pretty; 

%%
figure('Position',[289 444 1085 400]); hold on;
subplot(1,2,1); hold on; 
linetype = {'-',':'};
for i = [1 2]
    cx = cell2mat(compinfo(i,:)'); 
    cxLH = bootstrap_LHM(@nanmean,cx); 
    cxLH = (cxLH - min(cxLH(:,3)))./range(cxLH(:,3)); 
    plot((1:size(cxLH,1))'*0.02,nanmean(cell2mat(onsetalign(i,:)),2),linetype{i},'Color',clrs(1,:)); 
    patchwithnan((1:size(cxLH,1))*0.02,cxLH(:,1),cxLH(:,2),clrs(4,:)); 
    plot((1:size(cxLH,1))'*0.02,cxLH(:,3),linetype{i},'Color',clrs(4,:)); 
end
xlabel('ms')

subplot(1,2,2); hold on; 
for j = 1:4 % Task
    for i = [1 2] % Participant
        if ~isempty(info_reduction{i}{j})
            plot((j+0.2*i)*[1 1],info_reduction{i}{j}([1 3]),'Color',clrs(j,:)); 
            plot(j+.2*i,info_reduction{i}{j}(2),mrks{i},'Color',clrs(j,:)); 
        end
 
    end
end
xlim([1 5]); ylim([0 100]); %plot([1 5],[100 100],'k:'); 
ylabel('% var reduction reg1 -> reg2')
make_pretty; 
set(gca,'XTick',1.5:1:4.5); xticklabels(condnames)