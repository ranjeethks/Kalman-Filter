%Ranjeeth KS, University of Calgary
function [pos_r_plh r_bias]= update_state(method,h_constraint,dx,pos_r,x0,r_bias,pos_r_plh,vel_r)

if(method == 'carti')

        correction = [pos_r r_bias];   
        pos_r(1) = correction(1) + dx(1);
        pos_r(2) = correction(2) + dx(2);
        pos_r(3) = correction(3) + dx(3);
        r_bias   = correction(4) + dx(4);
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
        
end

        
       