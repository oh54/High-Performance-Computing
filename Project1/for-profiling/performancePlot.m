%% Load and plot data

% filename
importData

ii = find(submitmkn(:,1) > 0, 1);

submitmkn(1:ii-1,1:2) = submitmkn(1:ii-1,2:3);
submitmkn(:,3) = [];

p1 = plot(log2(submitmkn(:,1)),submitmkn(:,2),'*-');
set(gca,'xticklabel',[2^2, 2^4, 2^6, 2^8, 2^10/(2^10), 2^12/(2^10), 2^14/(2^10), 2^16/(2^10), 2^18/(2^10)])
ylim([1440 1600])
hold on
plot([log2(32) log2(32)],[-100, 1e5],'k--')
plot([log2(256) log2(256)],[-100, 1e5],'k--')
plot([log2(256*2^10) log2(256*2^10)],[-100, 1e5],'k--')