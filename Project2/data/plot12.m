%% plot for exercise 1 and 2 proejct 2
% loads time, N, d, k
clear all

addpath('/zhome/f7/e/77392/02614_HPC/Week1/Assignment1/High-Performance-Computing/Project2/data/');

filename = 'change_d.txt';

importData

changed_jac = changed(1:end/2,:);
changed_gs = changed(end/2+1:end,:);
%%
semilogx(changed_jac(:,3),changed_jac(:,4)*1e-3)
hold on
semilogx(changed_gs(:,3),changed_gs(:,4)*1e-3)
axis tight
hl = legend('Jacobi','Gauss-Seidel','Location','best');
set(hl,'FontSize',14)
set(gca,'FontSize',14)
set(gca,'Xtick',[1e-10 1e-8 1e-6 1e-4 1e-2 1e0])
xlabel('Threshold','FontSize',18)
ylabel('Iterations [thousands]','FontSize',18)