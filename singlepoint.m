%Ranjeeth KS, University of Calgary
function  [run]=singlepoint(InpData,settings)


%%
%Initializing task related variables

epochs=settings.epochs;  %All the epochs are considered
no_iterationsinLS= settings.no_iterationsinLS;
task4check = settings.task4check; 
method = settings.method;
h_constraint = settings.h_constraint;
num_states = settings.num_states; %number of states
alpha = settings.alpha;
beta = settings.beta;
weightedDOP = settings.weightedDOP;

%%
%True position of receiver
r_lat_true = 51.0 + 04/60 + 47.83343/3600;
r_long_true =-(114.0 + 08/60 + 1.85354/3600) ;
r_alt_true = 1118.792;
aPriori_var_factor = settings.sigma02;

%%
 gTestFail_count = 0;
 lTestFail_count = 0;
 
%k'th epoch
%sigmaarr = [6 1];
%for weightedDOP =0:1    
%aPriori_var_factor = sigmaarr(weightedDOP+1);
 for k=1:epochs
     
    number_of_sats = InpData.numsats(k); 
    current_comb=[1:number_of_sats];
    PRNsinepoch = InpData.satdata(k).prnID;
     for n=1:number_of_sats
                 current_prnInd = current_comb(n);                 
                 pos_s(n,:) =  InpData.satdata(k).sat_data(n,5:7);                 
                 z(n)=InpData.satdata(k).sat_data(n,1);               
                 true_el(n) = InpData.satdata(k).sat_data(n,4)*pi/180;
                 Qr_temp(n,n)=1/((sin(true_el(n))*sin(true_el(n)))^weightedDOP);
     end
     
 %repeat until blunder, if present, is detected and rejected for epoch k
 blunderYes = 1;  
 blunder_iter = 1;
 use_aPost_var_factor = 0;
 sigma02 = aPriori_var_factor;
 PRNIDarray = PRNsinepoch;
 elevationarray = true_el;
 
 while(blunderYes~=0)
 
     % Initial guess of receiver state vector 'x' in each epoch
     r_lat =51;%.079953730555555555555555555556;%40;
     r_long =-114;%.13384820555555555555555555556;
     r_alt = 1000;
     r_bias = 0;

     r_lat = r_lat*pi/180;
     r_long = r_long*pi/180;

     %Initial guess of  states and uncertainity
     x0 = [1118.7]; 
     P0 = [0.1] ;
 
    if(~h_constraint)
    pos_r_plh = [r_lat r_long r_alt];
    else
    pos_r_plh = [r_lat r_long x0];   
    end
    pos_r = plh2xyz(pos_r_plh(1), pos_r_plh(2), pos_r_plh(3)); 
    if(use_aPost_var_factor == 1)
       sigma02 = sigma02_hat;
    end
        
    %Begin iteration
     for iter=1:no_iterationsinLS

           %get dz, H, and Qr matrix..      
           [dz H]=compute_dzandH(method,z,pos_s,pos_r,r_bias,pos_r_plh,length(z));  
           Qr = Qr_temp;  
           %append height as an observation if using height constraint
           if(h_constraint)
               dx0 = x0-pos_r_plh(3);
               dz = [dz dx0];
               H = [H;
                   0 0 1 0];
               Qr = [Qr (zeros(1,length(Qr(:,1))))';
                     zeros(1,length(Qr(1,:))) P0/sigma02];
           end                
        R=sigma02*Qr;        
        %Least square computation
        dx = computeLS(h_constraint,dz', H, R);       
        %co-factor matrix of states       
        Qp = computeQp(h_constraint,H,Qr);

        %Position correction
        [pos_r_plh r_bias]= update_state(method,h_constraint,dx,pos_r,x0,r_bias,pos_r_plh);
         pos_r = plh2xyz(pos_r_plh(1), pos_r_plh(2), pos_r_plh(3)); 
         
         store_dz = dz;
         store_H = H;
         store_R = R;
         store_Qr = Qr;
         store_Qp = Qp;
         
         clear('dz','H','R','Qr','Qp','P');
     end
 
        %Compute error in curvilinear frame and translate to error in local co-ordinate frame  
         dp=(r_lat_true*pi/180-pos_r_plh(1));
         dl=(r_long_true*pi/180-pos_r_plh(2));
         dh=r_alt_true-pos_r_plh(3);

         N=calcn(pos_r_plh(1));
         M=calcm(pos_r_plh(1));

         de = dl*(N+pos_r_plh(3))*cos(pos_r_plh(1));
         dn = dp*(M+pos_r_plh(3));
         dv = dh;
 
         %Tranlate Qp matrix in ECEF to ENU
         if(method == 'carti')
             store_Qp = TransDOPtoENU(pos_r_plh(1),pos_r_plh(2),store_Qp);
         end
  
         store_P=sigma02*store_Qp;
         %compute residue 
         [r]=compute_dzandH(method,z,pos_s,pos_r,r_bias,pos_r_plh,length(z)); 
             if(h_constraint)
                 dx0 = x0-pos_r_plh(3);
                 r = [r dx0];               
             end

         %residue is a column vector
         r=r';
         Cr = store_R - store_H*inv(store_H'*inv(store_R)*store_H)*store_H';

         run(blunder_iter).residual(k).res=r';
         run(blunder_iter).PRNIDarray(k).prn=PRNIDarray;
         run(blunder_iter).elevationarray(k).ele=elevationarray;
 
         %a-posteriori variance factor
         DOF = length(z)-num_states; %degree of freedom
         sigma02_hat = r'*inv(store_Qr)*r/DOF;
 
                 %%%%%%%%%%%%%%% GLOBAL TEST %%%%%%%%%%%%%%% 
                 %set thresholds
                 X2_1_minus_alpha_by_2 = chi2inv(alpha/2,DOF);
                 X2_alpha_by_2 = chi2inv(1-alpha/2,DOF);                 
                 zetaG = r'*inv(store_R)*r;
                 run(blunder_iter).epoch(k).zetaG=zetaG;
                 %test code
                 % global_test = 'pass';
                 % blunderYes = 0;
                 if(X2_alpha_by_2> zetaG & zetaG > X2_1_minus_alpha_by_2)
                     global_test = 'pass';
                     run(blunder_iter).GTestFailepochs(k)=0;
                     run(blunder_iter).num_blunders(k) = 0;
                     blunderYes = 0;
                 else
                     global_test = 'fail';
                     run(blunder_iter).GTestFailepochs(k)=1;
                     gTestFail_count = gTestFail_count + 1;
                 end
 
                 if(global_test == 'fail')

                     %recompute R based on a-posteriori estimate of co-variance factor
                     store_R = sigma02_hat*store_Qr;
                     
                     %%%%%%%%%%%%%%% LOCAL TEST %%%%%%%%%%%%%%% 
                     Cr = store_R - store_H*inv(store_H'*inv(store_R)*store_H)*store_H';
                     %set threshold
                     N_1_minus_alpha_by_2 = norminv(1-alpha/2);
                     
                     for n=1:length(z)
                     zetaL(n)=abs(r(n)/sqrt(Cr(n,n)));
                     end
                     
                     lTestFailInd=find(zetaL > N_1_minus_alpha_by_2);
                     if(lTestFailInd)
                         lTestFail_count = lTestFail_count + 1;
                         blunderYes = 1;
                         use_aPost_var_factor = 1;
                     else
                         blunderYes = 0;
                     end

                     run(blunder_iter).num_blunders(k) = length(lTestFailInd);
                     run(blunder_iter).blunder_sat_index(k).ind = lTestFailInd;
                     
                         if(run(blunder_iter).num_blunders(k)>1)
                             [max_zetal max_ind ]=max(zetaL);
                             sat_ind_to_reject = max_ind;
                             run(blunder_iter).Sat_to_reject(k)=PRNIDarray(sat_ind_to_reject); 
                         elseif(run(blunder_iter).num_blunders(k) == 1)
                             sat_ind_to_reject = run(blunder_iter).blunder_sat_index(k).ind;
                             run(blunder_iter).Sat_to_reject(k)=PRNIDarray(sat_ind_to_reject); 

                         end
                         
                 end
                 
                 
                 %Internal and external reliability
                 %compute delta0
                 delta0 = norminv(1-alpha/2,0,1) + norminv(1-beta,0,1);
                 %compute int and ext reliability before rejecting any blunder
                 if(blunder_iter == 1)
                     for n=1:length(z)
                     MDB(k).sat(n) = delta0*store_R(n,n)/sqrt(Cr(n,n));
                     ext_rel = (inv(store_H'*inv(store_R)*store_H))*(store_H'*inv(store_R)*MDB(k).sat(n));
                     dx_ext_rel(k).Nsats(n,1)=sqrt(ext_rel(1).^2 + ext_rel(2).^2);
                     dx_ext_rel(k).Nsats(n,2)=ext_rel(3);
                     end
                 end
                 
        %%STORE OUTPUTS %%
          run(blunder_iter).dENU(k,1:3)=[de dn dv];
          run(blunder_iter).stdENU(k,1:4)= [sqrt(store_P(1,1)) sqrt(store_P(2,2)) sqrt(store_P(3,3)) sqrt(store_P(4,4))];
          run(blunder_iter).co_var(k,1:6)=[store_P(1,2) store_P(1,3) store_P(1,4) store_P(2,3) store_P(2,4) store_P(3,4)];
          run(blunder_iter).DOP(k,:)=computeDOP(store_Qp);
          run(blunder_iter).num_sats(k)=length(z);
 
          %plot error ellipse
          if(blunderYes==0 & (k == 1 | k==epochs))
              %plot 2D error ellipse
                xc = 0; yc = 0; % Center of ellip 
                sigma2E = store_P(1,1);
                sigma2N = store_P(2,2);
                sigmaNE = store_P(2,1);
                a = sqrt(0.5*(sigma2E + sigma2N)+sqrt(0.25*(sigma2E - sigma2N).^2+(sigmaNE).^2));   % Major and minor radius 
                b = sqrt(0.5*(sigma2E + sigma2N)-sqrt(0.25*(sigma2E - sigma2N).^2+(sigmaNE).^2));
                theta = 0.5*atan2(2*(sigmaNE),sigma2N-sigma2E);       % Orientation (degree) 
                phi = pi/2 - theta;
                num = 100;      % #points -> smoothness 
                angle = linspace(0,2*pi,num);
                p(1,:) = a*cos(angle);
                p(2,:) = b*sin(angle);
                Q = [cos(phi) -sin(phi)
                    sin(phi)  cos(phi)];
                p = Q*p;
                figure(1)
                plot(p(1,:),p(2,:),'g','LineWidth',2); hold on;
                plot(xc,yc,'r+');
                hold on;
          end
    
         %reject the observation if blunder is present
         if(global_test == 'fail' & blunderYes == 1)
             temp_z = z;
             temp_pos_s = pos_s;
             temp_Qr = Qr_temp;
             temp_PRNIDarray = PRNIDarray;
             clear('z','store_dz','store_H','store_R','store_Qr','store_Qp','store_P','pos_s','r','Cr','zetaL','lTestFailInd','Qr_temp','PRNIDarray');
             z = [temp_z(1:sat_ind_to_reject-1) temp_z(sat_ind_to_reject+1:end)];
             pos_s = [temp_pos_s(1:sat_ind_to_reject-1,:); temp_pos_s(sat_ind_to_reject+1:end,:)];
             Qr_temp = [temp_Qr(1:sat_ind_to_reject-1,1:sat_ind_to_reject-1) temp_Qr(1:sat_ind_to_reject-1,sat_ind_to_reject+1:end);
                        temp_Qr(sat_ind_to_reject+1:end,1:sat_ind_to_reject-1) temp_Qr(sat_ind_to_reject+1:end,sat_ind_to_reject+1:end)];
             blunder_iter = blunder_iter+1;
             PRNIDarray = [temp_PRNIDarray(1:sat_ind_to_reject-1) temp_PRNIDarray(sat_ind_to_reject+1:end)];
             clear('temp_z','temp_pos_s','temp_Qr');
         end 
 end
 clear('z','store_dz','store_H','store_R','store_Qr','store_Qp','store_P','pos_s','r','Cr','zetaL','lTestFailInd','lTestFailInd','Qr_temp');
 
end
 %%
% weighing(weightedDOP+1).dENU = run.dENU;
% weighing(weightedDOP+1).stdENU = run.stdENU;
% end
stop = 1;
 
 if(gTestFail_count/epochs > 0.03 & gTestFail_count/epochs < 0.07)
     disp('a priori variance factor is tuned');
     gTestFail_count*100/epochs
 else
    disp('tune a priori variance factor');
     gTestFail_count*100/epochs
 end
 
 [residue elevation availablePRN MDBout dx_ext_relH dx_ext_relV]=arrangeResidue(run,epochs,MDB,dx_ext_rel);
  
  %%
  
  %%
 