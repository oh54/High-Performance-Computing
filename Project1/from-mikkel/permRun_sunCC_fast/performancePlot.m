%% Load and plot data

names = {'nmk','nkm','kmn','knm','mkn','mnk'}

for i = 1:length(names)
	filename = ['/zhome/f7/e/77392/02614_HPC/Week1/Assignment1/High-Performance-Computing/Project1/from-mikkel/permRun_sunCC_fast/' names{i} '_o.txt'];

	importData

	ii = find(dat(:,1,i) > 0, 1);

	dat(1:ii-1,1:2,i) = dat(1:ii-1,2:3,i);
	%dat(:,3,i) = [];

end
%%
figure(1)
for j = 1:length(names)
p1 = plot(log2(dat(:,1,j)),dat(:,2,j),'*-');
set(gca,'xticklabel',[2^2, 2^4, 2^6, 2^8, 2^10/(2^10), 2^12/(2^10), 2^14/(2^10), 2^16/(2^10), 2^18/(2^10) 2^20/(2^10)])
ylim([140 1600])
hold on
end
hl = legend(names,'location','east');
plot([log2(32) log2(32)],[-100, 1e5],'k--')
plot([log2(256) log2(256)],[-100, 1e5],'k--')
plot([log2(256*2^10) log2(256*2^10)],[-100, 1e5],'k--')

ylabel('Performance [MFlop/s]','FontSize',18)
xlabel('Memory footprint [kB and MB]','FontSize',18)
set(gca,'FontSize',16)
set(hl,'FontSize',16)
