%Ranjeeth KS, University of Calgary
function [Out] = KalmanFilter(InpData, settings)

%Initializing task related variables

epochs=settings.epochs;  %All the epochs are considered
no_iterationsinLS= settings.no_iterationsinLS;
method = settings.method;
h_constraint = 0; %h_constraint is not used
num_states = settings.num_states; %number of states
weightedDOP = settings.weightedDOP;
states = settings.states;
measurement = settings.measurement;

c = 299792458;
h0 = 1e-21;
h2 = 1e-20;
qb = c*c*h0/2;
qd = 2*pi*pi*c*c*h2;

%Pos_state_spectal_density: All options are listed in this array.
PosSD = [(1).^2 (0.1).^2 (0.01).^2 (0.0).^2]; 

%Vel_state_spectal_density: All options are listed in this array.
VelSD = [(1).^2 (0.01).^2]; 

%number of states
number_of_states = num_states;

% Intial guess of position and bias states
r_lat =51;
r_long =-114;
r_alt = 1000;
r_bias = 0;
 
r_drift = 0;
initialPdrift = 100; %variance of drift

r_ve = 0;
initialPve = 400; % assumes nominal speed of about 80km/hr

r_vn = 0;
initialPvn = 400; % assumes nominal speed of about 80km/hr

r_vu = 0;
initialPvu = 10; %kept minimal, assumes least vertical motion

vel_r = [r_ve r_vn r_vu];
clock = [r_bias r_drift];
 
%intial guess of std. deviation of velocities
stdENU_vel(1,:)=[sqrt(initialPve) sqrt(initialPvn) sqrt(initialPvu)];

aPriori_var_factorPR = settings.sigma02;
aPriori_var_factorPRR = settings.sigma02prr;
 
%%
%True position of receiver
r_lat_true = 51.0 + 04/60 + 47.83343/3600;
r_long_true =-(114.0 + 08/60 + 1.85354/3600) ;
r_alt_true = 1118.792;

%read true kinematic position
fid = fopen('D:\Study\ENGO 620\Lab2\Code\ReferenceSolution.txt');
%fid = fopen('D:\Study\ENGO 620\Lab2\Code\SimulatorTruth.txt');

TrueTraj = textscan(fid, '%f %f %f %f %f %f %f','HeaderLines');
fclose(fid);
[row column]=size(TrueTraj{1,1});

%copy true kinematic position
k=1;
    for i = 1:row   
        GPS_time(k)=TrueTraj{1,1}(i);
        KinematicLat(k)=TrueTraj{1,2}(i);
        KinematicLong(k)=TrueTraj{1,3}(i);
        KinematicAlt(k)=TrueTraj{1,4}(i);
        KinematicVn(k)=TrueTraj{1,5}(i);
        KinematicVe(k)=TrueTraj{1,6}(i);
        KinematicVu(k)=TrueTraj{1,7}(i);
        k=k+1;
    end
 
%%
% LS in first epoch
 for k=1:1
     
    number_of_sats = InpData.numsats(k); 
    current_comb=[1:number_of_sats];
    PRNsinepoch = InpData.satdata(k).prnID;
    n=0;
     for nn=1:number_of_sats
         if(InpData.satdata(k).sat_data(nn,4)>=9)
             n=n+1;
                 current_time = InpData.gpstime(k);
                               
                 pos_s(n,:) =  InpData.satdata(k).sat_data(nn,5:7);                 
                 z(n)=InpData.satdata(k).sat_data(nn,1);               
                 true_el(n) = InpData.satdata(k).sat_data(nn,4)*pi/180;
                 Qr(n,n)=1/((sin(true_el(n))*sin(true_el(n)))^weightedDOP);
                 
                 vel_s(n,:)=InpData.satdata(k).sat_data(nn,8:10); 
                 zv(n)=InpData.satdata(k).sat_data(nn,2); 
         end
     end
    
    temp_storeQr = Qr;
    sigma02 = aPriori_var_factorPR;
    sigma02PRR = aPriori_var_factorPRR;
    
    PRNIDarray = PRNsinepoch;
    elevationarray = true_el;
        
     if(settings.states == 'PVC')
         if(GPS_time(k)~=current_time)
             stop=1;
         end
         r_lat_true = KinematicLat(k);
         r_long_true = KinematicLong(k);
         r_alt_true = KinematicAlt(k);
         r_ve_true = KinematicVe(k);
         r_vn_true = KinematicVn(k);
         r_vu_true = KinematicVu(k);            
     end

     % Initial guess of receiver state vector 'x' 
     r_lat =51;
     r_long =-114;
     r_alt = 1000;
     r_bias = 0;
     
     r_lat = r_lat*pi/180;
     r_long = r_long*pi/180;

     %Initial guess of  states and uncertainity
     x0 = [1118.7]; 
     P0 = [0.1] ;
 
    pos_r_plh = [r_lat r_long r_alt];   

    pos_r = plh2xyz(pos_r_plh(1), pos_r_plh(2), pos_r_plh(3)); 
      
    %Begin iteration
     for iter=1:no_iterationsinLS

       %get dz, H, and Qr matrix.. 
       [dz H]=compute_dzandH(method,z,pos_s,pos_r,r_bias,pos_r_plh,length(z));  
       [dzv]=compute_dzvKF(method,zv,pos_s,vel_s,pos_r,pos_r_plh,vel_r,clock,length(zv),states);

       clear('Qr');
       Qr=temp_storeQr;
        
        Rv = sigma02PRR*Qr;   
        R=sigma02*Qr;      

        if(measurement == 'PRR')
            dz = [dz dzv];
            temp9 = 0*eye(length(z));
            R = [R temp9;
            temp9 Rv];
            Qr = [Qr temp9;
            temp9 Qr];
            clear('temp9');
            temp8 = (zeros(length(z),4));
            H = [H temp8;
            temp8 H];
            clear('temp8');
        end

        %Least square computation
        dx = computeLS(h_constraint,dz', H, R);          
       
        %Position correction
        [pos_r_plh r_bias]= update_state(method,h_constraint,dx,pos_r,x0,r_bias,pos_r_plh );
         pos_r = plh2xyz(pos_r_plh(1), pos_r_plh(2), pos_r_plh(3)); 

         if(measurement == 'PRR')
              vel_r(1) = vel_r(1) + dx(5);
              vel_r(2) = vel_r(2) + dx(6);
              vel_r(3) = vel_r(3) + dx(7);
              r_drift = r_drift   + dx(8);
         end

     end
            
     Qp = computeQp(h_constraint,H,Qr);  
     
     %Compute error in curvilinear frame and translate to error in local co-ordinate frame
     dp=r_lat_true*pi/180-pos_r_plh(1);
     dl=r_long_true*pi/180-pos_r_plh(2);
     dh=r_alt_true-pos_r_plh(3);

     N=calcn(pos_r_plh(1));
     M=calcm(pos_r_plh(1));

     de = dl*(N+pos_r_plh(3))*cos(pos_r_plh(1));
     dn = dp*(M+pos_r_plh(3));
     dv = dh;

     %Tranlate Qp matrix in ECEF to ENU
     if(method == 'carti')
         Qp = TransDOPtoENU(pos_r_plh(1),pos_r_plh(2),Qp);
     end
            
     if(measurement == 'PRR') 
         
         P_matrix = [sigma02*Qp(1:4,:); sigma02PRR*Qp(5:8,:)];
         dve = r_ve_true - vel_r(1);
         dvn = r_vn_true - vel_r(2);
         dvu = r_vu_true - vel_r(3);
         stdEv = sqrt(P_matrix(5,5));
         stdNv = sqrt(P_matrix(6,6));
         stdUv = sqrt(P_matrix(7,7));
         dENU_vel(k,:)=[dve dvn dvu];             
         stdENU_vel(k,:)=[stdEv stdNv stdUv];
         
     elseif(measurement == 'PRO')
         
         P_matrix = sigma02*Qp;
         
     end
         
     stdE = sqrt(P_matrix(1,1));
     stdN = sqrt(P_matrix(2,2));
     stdU = sqrt(P_matrix(3,3));

     dENU(k,:)=[de dn dv];
     stdENU(k,:)=[stdE stdN stdU];
                            
 end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% L S Initialization complete %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% K F BEGINS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 %Set number of states based on the task
 if(settings.states == 'PCO')
     number_of_states = 5; % 4 states for testing purpose only
     CheckforPSD = PosSD;
 elseif(settings.states == 'PVC' )
     number_of_states = 8;
     CheckforPSD = VelSD;  
 end
     
% store apriori info from LS
if(number_of_states == 4) 
    
    KFaposterioriX = ([pos_r_plh r_bias  ])'; 
    KFaposterioriP = P_matrix;  
    
elseif(number_of_states == 5)  
    
    KFaposterioriX = ([pos_r_plh r_bias r_drift])';      
    temp = [0 0 0 0];
    KFaposterioriP = [P_matrix temp'; temp   initialPdrift]; 
    
elseif(number_of_states == 8 & measurement == 'PRO')     
    
    KFaposterioriX = ([pos_r_plh r_bias vel_r r_drift])'; 
    temp4 = zeros(4,4); 
    temp5 = [initialPve     0       0       0;
                0   initialPvn      0       0;
                0           0   initialPvu  0;
                0           0       0 initialPdrift];

    KFaposterioriP = [P_matrix temp4;
                      temp4    temp5];
elseif(number_of_states == 8 & measurement == 'PRR') 
    
    KFaposterioriX = ([pos_r_plh r_bias vel_r r_drift])'; 
    KFaposterioriP = P_matrix;
    
end

    % run the KF for all epochs, for different spectral density            
    for SD_iter = 1:length(CheckforPSD)
     
        if(settings.states == 'PCO')
            % if only 5 states are considered, processnoise is included
            % only for position, bias and drift states; velocity states are
            % not considered
            spectral_density_vector = [CheckforPSD(SD_iter) CheckforPSD(SD_iter) CheckforPSD(SD_iter) qb 0 0 0 qd]; % E N U bias drift
        elseif(settings.states == 'PVC')
            spectral_density_vector = [0 0 0 qb CheckforPSD(SD_iter) CheckforPSD(SD_iter) CheckforPSD(SD_iter) qd]; % E N U bias Ve Vu Vu drift
        end  
     
        current_time = InpData.gpstime(1);        
        estimatedposition(SD_iter).plh(1,:) = [pos_r_plh(1)*180/pi pos_r_plh(2)*180/pi pos_r_plh(3)];

            for k=2:epochs

            clear('dz','H','R','Qr','z','zv','dzv','Rv');
            previous_time = current_time;
            current_time = InpData.gpstime(k);
            dt = current_time-previous_time;

            number_of_sats = InpData.numsats(k); 
            current_comb=[1:number_of_sats];
            PRNsinepoch = InpData.satdata(k).prnID;        
            n=0;

             for nn=1:number_of_sats
                     n=n+1;
                     current_time = InpData.gpstime(k);

                     pos_s(n,:) =  InpData.satdata(k).sat_data(nn,5:7);                 
                     z(n)=InpData.satdata(k).sat_data(nn,1);               
                     true_el(n) = InpData.satdata(k).sat_data(nn,4)*pi/180;
                     Qr(n,n)=1/((sin(true_el(n))*sin(true_el(n)))^weightedDOP);

                     vel_s(n,:)=InpData.satdata(k).sat_data(nn,8:10); 
                     zv(n)=InpData.satdata(k).sat_data(nn,2); 
             end

            R=sigma02*Qr;
            Rv = sigma02PRR*Qr;

            if(settings.states == 'PVC')
                if(GPS_time(k)~=current_time)
                    stop=1;
                end
                r_lat_true = KinematicLat(k);
                r_long_true = KinematicLong(k);
                r_alt_true = KinematicAlt(k);
                r_ve_true = KinematicVe(k);
                r_vn_true = KinematicVn(k);
                r_vu_true = KinematicVu(k);            
            end

             TMatrix = GetTransitionMatrix(settings.states,dt,number_of_states);
             QMatrix = GetProcessNoiseMatrix(settings.states,dt,spectral_density_vector,number_of_states);

             if(settings.states == 'PCO')
                KFaprioriX = TMatrix*KFaposterioriX;
             elseif(settings.states == 'PVC')
                 %the below conversion in transition matrix is to predict latitude and longitude from North velocity and east velocity
                 %respectively. However this doesn't change Q matrix as all 'delta' position states are arranged in the order of their
                 %respective velocity states, and all are in the units of meters.
                N=calcn(KFaposterioriX(1));
                M=calcm(KFaposterioriX(1));    
                temp7 = TMatrix;
                temp7(1,6)=TMatrix(1,5)/(M+KFaposterioriX(3));
                temp7(1,5)= 0;
                temp7(2,5)=TMatrix(2,6)/((N+KFaposterioriX(3))*cos(KFaposterioriX(1)));
                temp7(2,6)= 0;
                KFaprioriX = temp7*KFaposterioriX;
             end

             KFaprioriP = TMatrix*KFaposterioriP*TMatrix' + QMatrix;         
             pos_r = plh2xyz(KFaprioriX(1), KFaprioriX(2), KFaprioriX(3)); 

             if(number_of_states == 5 )
                 pos_r_plh = [KFaprioriX(1) KFaprioriX(2) KFaprioriX(3)];
                 clock = [KFaprioriX(4) KFaprioriX(5)];
             elseif(number_of_states == 4)
                 pos_r_plh = [KFaprioriX(1) KFaprioriX(2) KFaprioriX(3)];
                 clock = [KFaprioriX(4) 0];
             elseif( number_of_states == 8)
                 pos_r_plh = [KFaprioriX(1) KFaprioriX(2) KFaprioriX(3)];
                 clock = [KFaprioriX(4) KFaprioriX(8)];
                 vel_r = [KFaprioriX(5) KFaprioriX(6) KFaprioriX(7)];
             end

             [dz H]=compute_dzandHKF(method,z,pos_s,pos_r,clock,pos_r_plh,length(z),states);
             [dzv]=compute_dzvKF(method,zv,pos_s,vel_s,pos_r,pos_r_plh,vel_r,clock,length(zv),states);

             if(measurement == 'PRR')
                 dz = [dz dzv];
                 temp9 = 0*eye(length(z));
                 R = [R temp9;
                     temp9 Rv];
                 clear('temp9');
             end

             if(number_of_states == 5)
                 temp3 = (zeros(1,length(z)))';
                 H = [H temp3];
                 clear('temp3');
             elseif(number_of_states == 4)
                 H = H;
             elseif( number_of_states == 8 & measurement == 'PRO')
                 temp6 = (zeros(length(z),4));
                 H = [H temp6];
                 clear('temp6');
             elseif( number_of_states == 8 & measurement == 'PRR')
                 temp8 = (zeros(length(z),4));
                 H = [H temp8;
                     temp8 H];
                 clear('temp8');
             end

             KMatrix = KFaprioriP*H'*(inv(H*KFaprioriP*H'+R));         
             dx = KMatrix*dz';        
             [pos_r_plh r_bias vel_r r_drift]= update_stateKF(method,h_constraint,dx,pos_r,x0,clock,pos_r_plh, states,number_of_states, vel_r);                  
             KFaposterioriP = (eye(length(KFaposterioriX))-KMatrix*H)*KFaprioriP;

             %computing error
             %Compute error in curvilinear frame and translate to error in local co-ordinate frame  
             dp=r_lat_true*pi/180-pos_r_plh(1);
             dl=r_long_true*pi/180-pos_r_plh(2);
             dh=r_alt_true-pos_r_plh(3);

             N=calcn(pos_r_plh(1));
             M=calcm(pos_r_plh(1));

             de = dl*(N+pos_r_plh(3))*cos(pos_r_plh(1));
             dn = dp*(M+pos_r_plh(3));
             dv = dh;

             %Tranlate Qp matrix in ECEF to ENU
             if(method == 'carti')
                 if(number_of_states == 4)
                 KFaposterioriP = TransDOPtoENUKF(pos_r_plh(1),pos_r_plh(2),KFaposterioriP,states);
                 end
             end

             stdE = sqrt(KFaposterioriP(1,1));
             stdN = sqrt(KFaposterioriP(2,2));
             stdU = sqrt(KFaposterioriP(3,3));

             dENU(k,:)=[de dn dv];             
             stdENU(k,:)=[stdE stdN stdU];

             if(number_of_states == 4)
                KFaposterioriX = ([pos_r_plh r_bias])';
             elseif(number_of_states == 5)   
                 KFaposterioriX = ([pos_r_plh r_bias r_drift])';
             elseif(number_of_states == 8)   
                KFaposterioriX = ([pos_r_plh r_bias vel_r r_drift])';
                 dve = r_ve_true - vel_r(1);
                 dvn = r_vn_true - vel_r(2);
                 dvu = r_vu_true - vel_r(3);
                 stdEv = sqrt(KFaposterioriP(5,5));
                 stdNv = sqrt(KFaposterioriP(6,6));
                 stdUv = sqrt(KFaposterioriP(7,7));
                 dENU_vel(k,:)=[dve dvn dvu];             
                 stdENU_vel(k,:)=[stdEv stdNv stdUv];
                 estimatedvelocity(SD_iter).venu(k,:) = vel_r;
             end

            estimatedposition(SD_iter).plh(k,:) = [pos_r_plh(1)*180/pi pos_r_plh(2)*180/pi pos_r_plh(3)];
                
        end
        run(SD_iter).dENU = dENU;
        run(SD_iter).stdENU = stdENU;

        if(settings.states == 'PVC')
            run(SD_iter).dENU_vel = dENU_vel;
            run(SD_iter).stdENU_vel = stdENU_vel;
        end
        
    end
 
     if(settings.states == 'PVC')
         runKFKinematic = run;
         Out = runKFKinematic;
     elseif(settings.states == 'PCO')
         runKFStatic = run;
         Out = runKFStatic;
     end
 %%
 close all;
 duration = 1:650;
 
 
 estE_in_meters = (estimatedposition(1).plh(duration,2) - estimatedposition(1).plh(1,2))*((N+pos_r_plh(3))*cos(pos_r_plh(1)));
 estN_in_meters = (estimatedposition(1).plh(duration,1) - estimatedposition(1).plh(1,1))*(M+pos_r_plh(3));
 
 trueE_in_meters = (KinematicLong(duration) - KinematicLong(1))*((N+pos_r_plh(3))*cos(pos_r_plh(1)));
 trueN_in_meters = (KinematicLat(duration) - KinematicLat(1))*(M+pos_r_plh(3));
 
plot(trueE_in_meters,trueN_in_meters,'b');
hold on;
plot(estE_in_meters,estN_in_meters,'r.');
 
 
 %%
figure(1);
t = 0:length(runKFKinematic(1).dENU(:,1))-1;

[ax,h1,h2]=plotyy(t,runKFKinematic(1).dENU(:,1),t,runKFKinematic(2).dENU(:,1));

xlabel('time epochs [arbitrary unit]','fontweight','bold','fontsize',20);

axes(ax(1)); ylabel('East error [m]','fontweight','bold','fontsize',20);
axes(ax(2)); ylabel('East error [m]','fontweight','bold','fontsize',20);
set(ax(1),'YLim',[0.5 2.5]);
set(ax(2),'YLim',[-40 40]);

title('East position error: For two different noise spectral density values','fontweight','bold','fontsize',20);

axes(ax(1)); legend('NSD = 1 m/s^2/\surdHz/axis','location','southwest');
axes(ax(2)); legend('NSD = 0.01 m/s^2/\surdHz/axis');

axes(ax(1)); set(gca,'fontweight','bold','fontsize',20);
axes(ax(2)); set(gca,'fontweight','bold','fontsize',20);
grid on;

%%
figure(2);
t = 0:length(runKFKinematic(1).dENU(:,2))-1


[ax,h1,h2]=plotyy(t,runKFKinematic(1).dENU(:,2),t,runKFKinematic(2).dENU(:,2));


xlabel('time epochs [arbitrary unit]','fontweight','bold','fontsize',20);

axes(ax(1)); ylabel('North error [m]','fontweight','bold','fontsize',20);
axes(ax(2)); ylabel('North error [m]','fontweight','bold','fontsize',20);
set(ax(1),'YLim',[-5 -3]);
set(ax(2),'YLim',[-50 50]);

title('North position error: For two different noise spectral density values','fontweight','bold','fontsize',20);

axes(ax(1)); legend('NSD = 1 m/s^2/\surdHz/axis','location','southwest');
axes(ax(2)); legend('NSD = 0.01 m/s^2/\surdHz/axis');

axes(ax(1)); set(gca,'fontweight','bold','fontsize',20);
axes(ax(2)); set(gca,'fontweight','bold','fontsize',20);

grid on;
%%
figure(3);
t = 0:length(runKFKinematic(1).dENU(:,3))-1


[ax,h1,h2]=plotyy(t,runKFKinematic(1).dENU(:,3),t,runKFKinematic(2).dENU(:,3));


xlabel('time epochs [arbitrary unit]','fontweight','bold','fontsize',20);

axes(ax(1)); ylabel('Vertical error [m]','fontweight','bold','fontsize',20);
axes(ax(2)); ylabel('Vertical error [m]','fontweight','bold','fontsize',20);
set(ax(1),'YLim',[-5 -1]);
set(ax(2),'YLim',[-20 20]);

title('Vertical position error: For two different noise spectral density values','fontweight','bold','fontsize',20);

axes(ax(1)); legend('NSD = 1 m/s^2/\surdHz/axis','location','southwest');
axes(ax(2)); legend('NSD = 0.01 m/s^2/\surdHz/axis');

axes(ax(1)); set(gca,'fontweight','bold','fontsize',20);
axes(ax(2)); set(gca,'fontweight','bold','fontsize',20);

grid on;
%%
%subplot(2,1,2)
figure(4);
t = 0:length(runKFKinematic(1).dENU_vel(:,1))-1

[ax,h1,h2]=plotyy(t,runKFKinematic(1).dENU_vel(:,1),t,runKFKinematic(2).dENU_vel(:,1));

xlabel('time epochs [arbitrary unit]','fontweight','bold','fontsize',20);

axes(ax(1)); ylabel('East error [m/s]','fontweight','bold','fontsize',20);
axes(ax(2)); ylabel('East error [m/s]','fontweight','bold','fontsize',20);
set(ax(1),'YLim',[-1.5 1.5]);
set(ax(2),'YLim',[-15 15]);

title('East velocity error: For two different noise spectral density values','fontweight','bold','fontsize',20);

axes(ax(1)); legend('NSD = 1 m/s^2/\surdHz/axis','location','southwest');
axes(ax(2)); legend('NSD = 0.01 m/s^2/\surdHz/axis');

axes(ax(1)); set(gca,'fontweight','bold','fontsize',20);
axes(ax(2)); set(gca,'fontweight','bold','fontsize',20);
grid on;

%%
figure(5);
t = 0:length(runKFKinematic(1).dENU_vel(:,2))-1


[ax,h1,h2]=plotyy(t,runKFKinematic(1).dENU_vel(:,2),t,runKFKinematic(2).dENU_vel(:,2));


xlabel('time epochs [arbitrary unit]','fontweight','bold','fontsize',20);

axes(ax(1)); ylabel('North error [m/s]','fontweight','bold','fontsize',20);
axes(ax(2)); ylabel('North error [m/s]','fontweight','bold','fontsize',20);
set(ax(1),'YLim',[-2 2]);
set(ax(2),'YLim',[-20 20]);

title('North velocity error: For two different noise spectral density values','fontweight','bold','fontsize',20);

axes(ax(1)); legend('NSD = 1 m/s^2/\surdHz/axis','location','southwest');
axes(ax(2)); legend('NSD = 0.01 m/s^2/\surdHz/axis');

axes(ax(1)); set(gca,'fontweight','bold','fontsize',20);
axes(ax(2)); set(gca,'fontweight','bold','fontsize',20);

grid on;
%%
figure(6);
t = 0:length(runKFKinematic(1).dENU_vel(:,3))-1


[ax,h1,h2]=plotyy(t,runKFKinematic(1).dENU_vel(:,3),t,runKFKinematic(2).dENU_vel(:,3));


xlabel('time epochs [arbitrary unit]','fontweight','bold','fontsize',20);

axes(ax(1)); ylabel('Vertical error [m/s]','fontweight','bold','fontsize',20);
axes(ax(2)); ylabel('Vertical error [m/s]','fontweight','bold','fontsize',20);
set(ax(1),'YLim',[-0.5 0.8]);
set(ax(2),'YLim',[-1.5 2.5]);

title('Vertical velocity error: For two different noise spectral density values','fontweight','bold','fontsize',20);

axes(ax(1)); legend('NSD = 1 m/s^2/\surdHz/axis','location','southwest');
axes(ax(2)); legend('NSD = 0.01 m/s^2/\surdHz/axis');

axes(ax(1)); set(gca,'fontweight','bold','fontsize',20);
axes(ax(2)); set(gca,'fontweight','bold','fontsize',20);

grid on;
%%
 
   
     
     
 
