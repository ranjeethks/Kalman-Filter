%Ranjeeth KS, University of Calgary
function analysis(run)






 %%
  data_log(1,1)=mean(dENU(:,1));
    data_log(2,1)=mean(runKFStatic(1).dENU(:,1));
    data_log(3,1)=mean(runKFStatic(2).dENU(:,1));
    data_log(4,1)=mean(runKFStatic(3).dENU(:,1));
    data_log(5,1)=mean(runKFStatic(4).dENU(:,1));
    
    data_log(1,2)=std(dENU(:,1));
    data_log(2,2)=std(runKFStatic(1).dENU(:,1));
    data_log(3,2)=std(runKFStatic(2).dENU(:,1));
    data_log(4,2)=std(runKFStatic(3).dENU(:,1));
    data_log(5,2)=std(runKFStatic(4).dENU(:,1));
    %%
    data_log1(1,1)=mean(dENU(:,2));
    data_log1(2,1)=mean(runKFStatic(1).dENU(:,2));
    data_log1(3,1)=mean(runKFStatic(2).dENU(:,2));
    data_log1(4,1)=mean(runKFStatic(3).dENU(:,2));
    data_log1(5,1)=mean(runKFStatic(4).dENU(:,2));
    
    data_log1(1,2)=std(dENU(:,2));
    data_log1(2,2)=std(runKFStatic(1).dENU(:,2));
    data_log1(3,2)=std(runKFStatic(2).dENU(:,2));
    data_log1(4,2)=std(runKFStatic(3).dENU(:,2));
    data_log1(5,2)=std(runKFStatic(4).dENU(:,2));
     %%
     
    data_log2(1,1)=mean(dENU(:,3));
    data_log2(2,1)=mean(runKFStatic(1).dENU(:,3));
    data_log2(3,1)=mean(runKFStatic(2).dENU(:,3));
    data_log2(4,1)=mean(runKFStatic(3).dENU(:,3));
    data_log2(5,1)=mean(runKFStatic(4).dENU(:,3));
    
    data_log2(1,2)=std(dENU(:,3));
    data_log2(2,2)=std(runKFStatic(1).dENU(:,3));
    data_log2(3,2)=std(runKFStatic(2).dENU(:,3));
    data_log2(4,2)=std(runKFStatic(3).dENU(:,3));
    data_log2(5,2)=std(runKFStatic(4).dENU(:,3));
    %%
figure(1);

plot(dENU(:,1),'m','linewidth',3); hold on;
plot(runKFStatic(1).dENU(:,1),'b','linewidth',3);hold on;
plot(runKFStatic(2).dENU(:,1),'r','linewidth',3);hold on;
plot(runKFStatic(3).dENU(:,1),'g','linewidth',3);hold on;
plot(runKFStatic(4).dENU(:,1),'k','linewidth',3);hold on;

xlabel('time epochs [arbitrary unit]','fontweight','bold','fontsize',20);
ylabel('East error [m]','fontweight','bold','fontsize',20);
title('East error: LS vs. KF for different process noise','fontweight','bold','fontsize',20);
legend('LS Solution','q = 1 m/s/\surdHz/axis','q = 0.1 m/s/\surdHz/axis','q = 0.01 m/s/\surdHz/axis','q = 0.0 m/s/\surdHz/axis');
set(gca,'fontweight','bold','fontsize',20);
grid on;

 figure(2);

plot(dENU(:,2),'m','linewidth',3); hold on;
plot(runKFStatic(1).dENU(:,2),'b','linewidth',3);hold on;
plot(runKFStatic(2).dENU(:,2),'r','linewidth',3);hold on;
plot(runKFStatic(3).dENU(:,2),'g','linewidth',3);hold on;
plot(runKFStatic(4).dENU(:,2),'k','linewidth',3);hold on;

xlabel('time epochs [arbitrary unit]','fontweight','bold','fontsize',20);
ylabel('North error [m]','fontweight','bold','fontsize',20);
title('North error: LS vs. KF for different process noise','fontweight','bold','fontsize',20);
legend('LS Solution','q = 1 m/s/\surdHz/axis','q = 0.1 m/s/\surdHz/axis','q = 0.01 m/s/\surdHz/axis','q = 0.0 m/s/\surdHz/axis');
set(gca,'fontweight','bold','fontsize',20);
grid on;

figure(3);

plot(dENU(:,3),'m','linewidth',3); hold on;
plot(runKFStatic(1).dENU(:,3),'b','linewidth',3);hold on;
plot(runKFStatic(2).dENU(:,3),'r','linewidth',3);hold on;
plot(runKFStatic(3).dENU(:,3),'g','linewidth',3);hold on;
plot(runKFStatic(4).dENU(:,3),'k','linewidth',3);hold on;

xlabel('time epochs [arbitrary unit]','fontweight','bold','fontsize',20);
ylabel('Up error [m]','fontweight','bold','fontsize',20);
title('Up error: LS vs. KF for different process noise','fontweight','bold','fontsize',20);
legend('LS Solution','q = 1 m/s/\surdHz/axis','q = 0.1 m/s/\surdHz/axis','q = 0.01 m/s/\surdHz/axis','q = 0.0 m/s/\surdHz/axis');
set(gca,'fontweight','bold','fontsize',20);
grid on;
%%