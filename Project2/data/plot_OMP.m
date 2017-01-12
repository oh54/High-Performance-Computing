clear all
close all
filename = 'part_o.txt';

importData_omp
% changed_jac = parto(1:end/2,:);
% changed_gs = c(end/2+1:end,:);

% speed up plot
f = 0.22;
figure(1)
plot(parto(:,end),parto(1,2)./parto(:,2),'*-')
hold on
plot(1./((1-f)+f./parto(:,end)),'k--')

hl = legend('Initial parallel',['Amdahl, f = ' num2str(f)],'Location','SouthEast');
set(gca,'FontSize',14)
set(hl,'FontSize',14)
xlabel('No. of processors','FontSize',18)
ylabel('Speed up','FontSize',18)



% efficiency plot
figure(2)
plot(parto(:,end),parto(1,2)./parto(:,2)./parto(:,end),'-*')
set(gca,'FontSize',14)
xlabel('No. of processors','FontSize',18)
ylabel('Efficiency','FontSize',18)

filename = 'part_o_omp2.txt';
importData_omp
hold on
plot(parto(:,end),parto(1,2)./parto(:,2)./parto(:,end),'-*')
hl = legend('omp','omp2');
set(hl,'FontSize',14)

figure(3)
f = 0.93;
plot(parto(:,end),parto(1,2)./parto(:,2),'*-')
hold on
plot(1./((1-f)+f./parto(:,end)),'k--')

hl = legend('Initial parallel',['Amdahl, f = ' num2str(f)],'Location','SouthEast');
set(gca,'FontSize',14)
set(hl,'FontSize',14)
xlabel('No. of processors','FontSize',18)
ylabel('Speed up','FontSize',18)
%%
N = 3;
filename = 'parn_o.txt';
importData_parn
changed_seq = parno(1:end/N,:);
changed_omp = parno(end/N+1:end/N*2,:);
changed_omp2 = parno(end/N*2+1:end,:);

figure(4)
plot(changed_seq(:,3),changed_seq(:,2)./changed_seq(:,4))
hold on
plot(changed_omp(:,3),changed_omp(:,2)./changed_omp(:,4))
plot(changed_omp2(:,3),changed_omp2(:,2)./changed_omp2(:,4))

