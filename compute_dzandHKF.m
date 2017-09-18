%Ranjeeth KS, University of Calgary
function [dz H]=compute_dzandHKF(method,z,pos_s,pos_r,clock,pos_r_plh,num_sats,states);  

r_bias = clock(1);
r_drift = clock(2);


    for s=1:num_sats
    relative_pos_vect = pos_s(s,:) - pos_r;

        if(method == 'carti') %single point
            p = sqrt(relative_pos_vect(1).^2+relative_pos_vect(2).^2+relative_pos_vect(3).^2);    

            H(s,:) = [-relative_pos_vect(1)/p -relative_pos_vect(2)/p -relative_pos_vect(3)/p 1 ];    

        elseif(method == 'curvi') %single difference
            x_ENU=ECEF2ENU(pos_r_plh(1),pos_r_plh(2),pos_r_plh(3),relative_pos_vect');
            [az el] = calcAzEl(x_ENU);

            p = sqrt(x_ENU(1).^2+x_ENU(2).^2+x_ENU(3).^2);

            H(s,:) = [-cos(el)*sin(az) -cos(el)*cos(az) -sin(el) 1 ];
        end

    z_est = p+r_bias;
    dz(s) = z(s)-z_est;
    end
