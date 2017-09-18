%Ranjeeth KS, University of Calgary
function [dzv]=compute_dzvKF(method,zv,pos_s,vel_s,pos_r,pos_r_plh,vel_r_ENU,clock,num_sats,states)

r_bias = clock(1);
r_drift = clock(2);


vel_r_ECEF = ENU2ECEF(pos_r_plh(1),pos_r_plh(2),pos_r_plh(3),vel_r_ENU');


    for s=1:num_sats
    relative_pos_vect = pos_s(s,:) - pos_r;
    row = sqrt(relative_pos_vect(1).^2+relative_pos_vect(2).^2+relative_pos_vect(3).^2);    
    
    rowdot = (pos_s(s,1)-pos_r(1))*(vel_s(s,1)-vel_r_ECEF(1))+(pos_s(s,2)-pos_r(2))*(vel_s(s,2)-vel_r_ECEF(2))+(pos_s(s,3)-pos_r(3))*(vel_s(s,3)-vel_r_ECEF(3));
  
    rowdot = rowdot/row;   

    z_est = rowdot+r_drift;
    dzv(s) = zv(s)-z_est;
    end
