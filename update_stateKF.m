%Ranjeeth KS, University of Calgary
function [pos_r_plh r_bias vel_r r_drift]= update_stateKF(method,h_constraint,dx,pos_r,x0,clock,pos_r_plh, states,number_of_states,vel_r)


r_bias = clock(1);
r_drift = clock(2);


    if(method == 'carti')

            correction = [pos_r r_bias];   
            pos_r(1) = correction(1) + dx(1);
            pos_r(2) = correction(2) + dx(2);
            pos_r(3) = correction(3) + dx(3);
            r_bias   = correction(4) + dx(4);
            
            if(number_of_states==5)
                r_drift = r_drift + dx(5);
            end
            %x,y,z to phi,lambda,height
            pos_r_plh = xyz2plh(pos_r);
            pos_r_plh = pos_r_plh';

    elseif(method == 'curvi')


     N=calcn(pos_r_plh(1));
     M=calcm(pos_r_plh(1));

     % translate error in local co-ordinate frame  
     dl = dx(1)/((N+pos_r_plh(3))*cos(pos_r_plh(1)));
     dp = dx(2)/(M+pos_r_plh(3));
     dv = dx(3);

            correction = [pos_r_plh r_bias];
            pos_r_plh(1) = correction(1) + dp;
            pos_r_plh(2) = correction(2) + dl;
            pos_r_plh(3) = correction(3) + dv;
            r_bias       = correction(4) + dx(4);
            if(number_of_states==5)
                r_drift = r_drift + dx(5);
            end
            
            if(number_of_states==8)
                vel_r(1) = vel_r(1) + dx(5);
                vel_r(2) = vel_r(2) + dx(6);
                vel_r(3) = vel_r(3) + dx(7);
                r_drift = r_drift   + dx(8);
            end
            
                

    end


        
       