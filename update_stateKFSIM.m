%Ranjeeth KS, University of Calgary
function [pos_r_plh r_bias vel_r r_drift thetaDot]= update_stateKFSIM(method,h_constraint,dx,pos_r,x0,clock,pos_r_plh, states,number_of_states,vel_r,thetaDot)


r_bias = clock(1);
r_drift = clock(2);




    if(method == 'curvi')


     N=calcn(pos_r_plh(1));
     M=calcm(pos_r_plh(1));

     % translate error in local co-ordinate frame  
     dl = dx(1)/((N+pos_r_plh(3))*cos(pos_r_plh(1)));
     dp = dx(2)/(M+pos_r_plh(3));
     dv = 0;

            correction = [pos_r_plh r_bias];
            pos_r_plh(1) = correction(1) + dp;
            pos_r_plh(2) = correction(2) + dl;
            pos_r_plh(3) = correction(3) + dv;
            r_bias       = correction(4) + dx(3);
            
            if(number_of_states==7)
                vel_r(1) = vel_r(1) + dx(4);
                vel_r(2) = vel_r(2) + dx(5);
                vel_r(3) = vel_r(3) + 0;
                r_drift = r_drift   + dx(6);
                thetaDot = thetaDot + dx(7);
                
            end
            
                

    end


        
       