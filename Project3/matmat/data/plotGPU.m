%% plot GPU
clear all
loadCPU
loadGPU

gpu1 = gpu15lib(1:8,:);


for i = 1:5
    gpu(:,:,i) = gpu15lib(9+16*(i-1):9+15+ 16*(i-1),:);
end

%%
close all
figure(1)
semilogy(log2(sqrt(cpuDGEMM(1:length(gpu1),1)/8/3*1024)),cpuDGEMM(1:length(gpu1),2),'*-')
hold on
semilogy(log2(sqrt(gpu1(:,1)/8/3*1024)),gpu1(:,2),'*-')
vec = [4 5 6 7 8];
set(gca,'xtick',vec)
set(gca,'xticklabel',2.^vec)
xlabel('Side of square matrix','FontSize',16)
ylabel('Performance [MFlops/s]','FontSize',16)
set(gca,'FontSize',14)
hl = legend('cblasDgemm','gpu1');
set(hl,'FontSize',14);

for j = 1:5;
figure(j+1)
semilogy(log2(sqrt(cpuDGEMM(:,1)/8/3*1024)),cpuDGEMM(:,2)*1e-3,'*-')
hold on
    for k = 1:j
        semilogy(log2(sqrt(gpu(:,1,k)/8/3*1024)),gpu(:,2,k)*1e-3,'*-')
        hold on
        legs{k+1} = ['gpu' num2str(k+1)];
    end
    legs{1} = 'cblasDgemm';
vec = 4:ceil(log2(gpu(end,1,1)));
set(gca,'xtick',vec)
set(gca,'xticklabel',2.^vec)
xlabel('Side of square matrix','FontSize',16)
ylabel('Performance [GFlops/s]','FontSize',16)
set(gca,'FontSize',14)
if j == 5
    legs{k+1} = 'cublasDgemm';
end
hl = legend(legs,'location','SouthEast');
set(hl,'FontSize',14);
grid on
end
