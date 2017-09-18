%Ranjeeth KS, University of Calgary
function [Out] = KalmanFilter(InpData, settings)

%Initializing task related variables

epochs=settings.epochs;  %All the epochs are considered
no_iterationsinLS= settings.no_iterationsinLS;
method = settings.method;
h_constraint = 0;
num_states = settings.num_states; %number of states
weightedDOP = settings.weightedDOP;
states = settings.states;
measurement = settings.measurement;

%spectral_density_vector = settings.spectral_density_vector;
c = 299792458;
h0 = 1e-21;
h2 = 1e-20;
qb = c*c*h0/2;
qd = 2*pi*pi*c*c*h2;

%number of states
number_of_states = num_states;

% Intial guess of states, these will be inputs to the LS, which in turn
% will generate a-prioti information to KF in first epoch
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

thetaDot = 0; % agnular rate of rotation
initialPtd = 1; % rad/s : intial uncertainity in angular velocity
 
%intial guess of std of velocities
stdENU_vel(1,:)=[sqrt(initialPve) sqrt(initialPvn) sqrt(initialPvu)];
stdthetaDot(1) = sqrt(initialPtd);

aPriori_var_factorPR = settings.sigma02;
aPriori_var_factorPRR = settings.sigma02prr;
 
%%
%True position of receiver
r_lat_true = 51.0 + 04/60 + 47.83343/3600;
r_long_true =-(114.0 + 08/60 + 1.85354/3600) ;
r_alt_true = 1118.792;

%read true kinematic position
fid = fopen('D:\Study\ENGO 620\Lab2\Code\SimulatorTruth.txt');
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
        %North and east velocities are swapped
        KinematicVe(k)=TrueTraj{1,5}(i);
        KinematicVn(k)=TrueTraj{1,6}(i);
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
     for iter=1:2

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
        
         %update velocity if LS Measurement model includes velocity and
         %drift states wlong with position and bias states. 
         %If only pseudo-range measurement is used, LS cannot update velocity
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
 
   
% store apriori info from LS
if(number_of_states == 7) 
    pos_r_plh = [pos_r_plh(1:2) r_alt_true]; %True altitude is considered
    vel_r = [vel_r(1:2) 0];            %Vertical velocity is assumed zero
    
    KFaposterioriX = ([pos_r_plh(1:2) r_bias vel_r(1:2) r_drift thetaDot])'; 
    KFaposterioriP = [P_matrix(1:2,1:2) P_matrix(1:2,4:6) P_matrix(1:2,8) zeros(2,1);
                      P_matrix(4:6,1:2) P_matrix(4:6,4:6) P_matrix(4:6,8) zeros(3,1);
                      P_matrix(8,1:2) P_matrix(8,4:6) P_matrix(8,8) zeros(1,1);
                      0     0       0       0       0       0   initialPtd];  
    
end

    % run the KF            
        SD_iter = 1;
            
        spectral_density_vector = [0 0 qb 40 40 qd 0.01]; % E N bias Ve Vu drift thetaDot
            
        current_time = InpData.gpstime(1);        
        estimatedposition(SD_iter).plh(1,:) = [pos_r_plh(1)*180/pi pos_r_plh(2)*180/pi 1000];
        
        for k=2:epochs
         
        clear('dz','H','R','Qr','z','zv','dzv','Rv');
        previous_time = current_time;
        current_time = InpData.gpstime(k);
        
        %Get time difference between the epochs
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
                disp('Signal time and reference time mismatch!');
                break;
            end
            r_lat_true = KinematicLat(k);
            r_long_true = KinematicLong(k);
            r_alt_true = KinematicAlt(k);
            r_ve_true = KinematicVe(k);
            r_vn_true = KinematicVn(k);
            r_vu_true = KinematicVu(k);            
        end
         
         if(number_of_states==7)
         % theta is tan inverse of -ve and vn, resolved in four quandrants    
         theta = atan2(-KFaposterioriX(4),KFaposterioriX(5));
         % speed is the magnitude of vector sum of ve and vn
         ConstSpeed = sqrt(KFaposterioriX(4).^2+KFaposterioriX(5).^2);
         
         %Get transition matrix
         TMatrix = GetTransitionMatrixSIM(settings.states,dt,number_of_states, theta, ConstSpeed);
         
         %Get process noise matrix
         QMatrix = GetProcessNoiseMatrixSIM(settings.states,dt,spectral_density_vector,number_of_states, theta, ConstSpeed);

             %the below conversion in transition matrix is to predict latitude and longitude from North velocity and east velocity
             %respectively. However this doesn't change Q matrix as all 'delta' position states are arranged in the order of their
             %respective velocity states, and all position states are in the units of meters.
            N=calcn(KFaposterioriX(1));
            M=calcm(KFaposterioriX(1));    
            temp7 = TMatrix;
            temp7(1,5)=TMatrix(1,4)/(M+r_alt_true);
            temp7(1,4)= 0;
            temp7(2,4)=TMatrix(2,5)/((N+r_alt_true)*cos(KFaposterioriX(1)));
            temp7(2,5)= 0;
            temp7(1,7)= -(ConstSpeed*sin(theta)*dt*dt/2)/(M+r_alt_true);
            temp7(2,7)= -(ConstSpeed*cos(theta)*dt*dt/2)/((N+r_alt_true)*cos(KFaposterioriX(1)));
            KFaprioriX = temp7*KFaposterioriX;
            
         KFaprioriP = TMatrix*KFaposterioriP*TMatrix' + QMatrix;   
         
         pos_r = plh2xyz(KFaprioriX(1), KFaprioriX(2), r_alt_true); % altitude = true_altitude
         pos_r_plh = (xyz2plh(pos_r))';
         clock = [KFaprioriX(3) KFaprioriX(6)];
         vel_r = [KFaprioriX(4) KFaprioriX(5) 0]; % north velocity = 0
         thetaDot = KFaprioriX(7);
         
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
          
         if( number_of_states == 7 & measurement == 'PRR')
             temp8 = (zeros(length(z),3));
             temp9 = H(1:length(z),1:2);
             temp10 = H(1:length(z),4);
             clear('H');
             H = [temp9 temp10 temp8 zeros(length(z),1);
                 temp8 temp9 temp10 zeros(length(z),1)];
            
             clear('temp8','temp9','temp10');
         end
         
         KMatrix = KFaprioriP*H'*(inv(H*KFaprioriP*H'+R));         
         dx = KMatrix*dz';      
                        
         [pos_r_plh r_bias vel_r r_drift thetaDot]= update_stateKFSIM(method,h_constraint,dx,pos_r,x0,clock,pos_r_plh, states,number_of_states, vel_r,thetaDot);                  
         KFaposterioriP = (eye(length(KFaposterioriX))-KMatrix*H)*KFaprioriP;
         
         theta_dot(k)=thetaDot;
         
         %computing error
         dp=r_lat_true*pi/180-pos_r_plh(1);
         dl=r_long_true*pi/180-pos_r_plh(2);

         N=calcn(pos_r_plh(1));
         M=calcm(pos_r_plh(1));

         de = dl*(N+pos_r_plh(3))*cos(pos_r_plh(1));
         dn = dp*(M+pos_r_plh(3));

         %Tranlate Qp matrix in ECEF to ENU
         if(method == 'carti')
             if(number_of_states == 4)
             KFaposterioriP = TransDOPtoENUKF(pos_r_plh(1),pos_r_plh(2),KFaposterioriP,states);
             end
         end
             
         stdE = sqrt(KFaposterioriP(1,1));
         stdN = sqrt(KFaposterioriP(2,2));
   
             if(number_of_states == 7)
                dENU(k,:)=[de dn 0];             
                stdENU(k,:)=[stdE stdN 0];               
                KFaposterioriX = ([pos_r_plh(1:2) r_bias vel_r(1:2) r_drift thetaDot])';
                dve = r_ve_true - vel_r(1);
                dvn = r_vn_true - vel_r(2);
                stdEv = sqrt(KFaposterioriP(5,5));
                stdNv = sqrt(KFaposterioriP(6,6));
                dENU_vel(k,:)=[dve dvn 0 ];             
                stdENU_vel(k,:)=[stdEv stdNv 0];
                estimatedvelocity(SD_iter).venu(k,:) = vel_r;
                estimatedposition(SD_iter).plh(k,:) = [pos_r_plh(1)*180/pi pos_r_plh(2)*180/pi pos_r_plh(3)];
             end
       
        end
        run(SD_iter).dENU = dENU;
        run(SD_iter).stdENU = stdENU;

        if(settings.states == 'PVC')
            run(SD_iter).dENU_vel = dENU_vel;
            run(SD_iter).stdENU_vel = stdENU_vel;
        end
 
     if(settings.states == 'PVC')
         runKFSimulation = run;
         Out = runKFSimulation;
     elseif(settings.states == 'PCO')
         runKFSimulation = run;
         Out = runKFSimulation;
     end
 %%
 
 clear('dz','H','R','Qr','z','zv','dzv','Rv','R','Qp');
 %%
 close all;
 duration = 1:650;
 
 
 estE_in_meters = (estimatedposition(1).plh(duration,2) - estimatedposition(1).plh(1,2))*((N+pos_r_plh(3))*cos(pos_r_plh(1)));
 estN_in_meters = (estimatedposition(1).plh(duration,1) - estimatedposition(1).plh(1,1))*(M+pos_r_plh(3));
 
 trueE_in_meters = (KinematicLong(duration) - KinematicLong(1))*((N+pos_r_plh(3))*cos(pos_r_plh(1)));
 trueN_in_meters = (KinematicLat(duration) - KinematicLat(1))*(M+pos_r_plh(3));
 
plot(trueE_in_meters,trueN_in_meters,'b.');
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

 
   
     
     
 
