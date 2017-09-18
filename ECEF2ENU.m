%Ranjeeth KS, University of Calgary
function x_ENU = ECEF2ENU(r_lat,r_long,r_alt,vector)

%
R_ecef2enu=[-sin(r_long)             cos(r_long)            0;
            -sin(r_lat)*cos(r_long) -sin(r_lat)*sin(r_long) cos(r_lat);
             cos(r_lat)*cos(r_long)  cos(r_lat)*sin(r_long) sin(r_lat)];

x_ENU = R_ecef2enu*vector;     