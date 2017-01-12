%% plot for exercise 1 and 2 proejct 2
% loads time, N, d, k
clear all

addpath('/zhome/f7/e/77392/02614_HPC/Week1/Assignment1/High-Performance-Computing/Project2/data/');

filename = 'seqd_o.txt';

importData

changed_jac = changed(1:end/2,:);
changed_gs = changed(end/2+1:end,:);
%%
figure(1)
semilogx(changed_jac(:,3),changed_jac(:,4)*1e-3,'-*')
hold on
semilogx(changed_gs(:,3),changed_gs(:,4)*1e-3,'*-')
semilogx(changed_jac(:,3),changed_jac(:,4)*1e-3/2,'k--')
axis tight
hl = legend('Jacobi','Gauss-Seidel','Jacobi divided by 2','Location','best');
set(hl,'FontSize',14)
set(gca,'FontSize',14)
set(gca,'Xtick',[1e-10 1e-8 1e-6 1e-4 1e-2 1e0])
xlabel('Threshold','FontSize',18)
ylabel('Iterations [thousands]','FontSize',18)


%% change n
clear all
filename = 'seqn_o_d01.txt';

importData
changed_jac = changed(1:end/2,:);
changed_gs = changed(end/2+1:end,:);


figure(2)
plot(log2(changed_jac(:,2).^2*3*8/1024),changed_jac(:,5).*changed_jac(:,2).^2./changed_jac(:,1)*1e-6,'*-');
hold on
plot(log2(changed_gs(:,2).^2*3*8/1024),changed_gs(:,5).*changed_gs(:,2).^2./changed_gs(:,1)*1e-6,'*-');
hold on
ylim([0 400])
plot([log2(32) log2(32)],[-100, 1e5],'k--')
plot([log2(256) log2(256)],[-100, 1e5],'k--')
plot([log2(25600) log2(25600)],[-100, 1e5],'k--')
set(gca,'xticklabel',[4 16 64 256 1024 4 16 64 256 1024])
xlabel('Memory footprint [kB and MB]','FontSize',25)
ylabel('Point update frequency [MHz]','FontSize',25)
hl = legend('Jacobi sequential','Gauss-Seidel sequential','Location','NorthEast');
set(hl,'FontSize',20)
set(gca,'FontSize',20)



figure(3)
loglog(changed_jac(:,2),changed_jac(:,4),'*-')
hold on
loglog(changed_gs(:,2),changed_gs(:,4),'*-')
ylim([0 6000])
%set(gca,'xticklabel',[4 16 64 256 1024 4 16 64 256 1024])

filename = 'seqn_o_d001.txt';
importData
changed_jac = changed(1:end/2,:);
changed_gs = changed(end/2+1:end,:);

loglog(changed_jac(:,2),changed_jac(:,4),'*-')
hold on
loglog(changed_gs(:,2),changed_gs(:,4),'*-')

xlabel('Side length of matrix','FontSize',18)
ylabel('Iterations to convergence','FontSize',18)
hl = legend('Jacobi seq, d = 0.01','Gauss-Seidel = d = 0.01','Jacobi seq, d = 0.001','Gauss-Seidel = d = 0.001','Location','SouthWest');
set(hl,'FontSize',16)
set(gca,'FontSize',16)