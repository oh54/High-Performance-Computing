%% compares different compiler options

addpath('/zhome/f7/e/77392/02614_HPC/Week1/Assignment1/High-Performance-Computing/Project1/from-mikkel')
runIDs = {'none','fast','xO5','xprefetch_level=3','unroll=8'};
legrunIDs = runIDs;
legrunIDs{4} = 'xprefetch\_level=3';
names = {'mkn','mnk','knm'};
for j = 1:length(names)
    for i = 1:length(runIDs)
        filename = ['/zhome/f7/e/77392/02614_HPC/Week1/Assignment1/High-Performance-Computing/Project1/from-mikkel/permRun_sunCC_' runIDs{i} '/' names{j} '_o.txt']; 
        importData
    end
    
    figure(j)
    for k = 1:length(runIDs)
        plot(log2(dat(:,1,k)),dat(:,2,k),'-*')
        hold on
    end
    hl = legend(legrunIDs,'location','best');
    set(gca,'xticklabel',[2^2, 2^4, 2^6, 2^8, 2^10/(2^10), 2^12/(2^10), 2^14/(2^10), 2^16/(2^10), 2^18/(2^10), 2^20/(2^10)])
    ylim([min(min(dat(:,2,:)))*0.8 max(max(dat(:,2,1:6)))*1.2])
    
    ylabel('Performance [MFlop/s]','FontSize',25)
    xlabel('Memory footprint [kB or MB]','FontSize',25)

    set(gca,'FontSize',20)
    set(hl,'FontSize',20)
end
