%% Load and plot data
clear dat
% filename
runID = 'xO5';
% runID = 'xpr';

names={'kmn','mkn','nmk','mnk','knm','nkm','lib','blk'};
addpath('/zhome/f7/e/77392/02614_HPC/Week1/Assignment1/High-Performance-Computing/Project1/from-mikkel')
for i = 1:length(names)
    filename = ['/zhome/f7/e/77392/02614_HPC/Week1/Assignment1/High-Performance-Computing/Project1/from-mikkel/permRun_sunCC_' runID '/' names{i} '_o.txt']; 
    importData
    
end
%%
close
for j = 1:6
    p1 = plot(log2(dat(:,1,j)),dat(:,2,j),'*-');
    set(gca,'xticklabel',[2^2, 2^4, 2^6, 2^8, 2^10/(2^10), 2^12/(2^10), 2^14/(2^10), 2^16/(2^10), 2^18/(2^10), 2^20/(2^10)])
    ylim([min(min(dat(:,2,:)))*0.8 max(max(dat(:,2,1:6)))*1.2])
    hold on
end
for k = 1:j
    legnames{k} = names{k};
end
hl = legend(legnames,'Location','east');
ylabel('Performance [MFlop/s]','FontSize',25)
xlabel('Memory footprint [kB or MB]','FontSize',25)

set(gca,'FontSize',20)
set(hl,'FontSize',20)

    plot([log2(32) log2(32)],[-100, 1e5],'k--')
    plot([log2(256) log2(256)],[-100, 1e5],'k--')
    plot([log2(25600) log2(25600)],[-100, 1e5],'k--')