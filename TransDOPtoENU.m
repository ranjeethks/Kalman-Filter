%Ranjeeth KS, University of Calgary
function Qp_ENU = TransDOPtoENU(r_lat,r_long,Qp_ECEF)

%Transformation matrix from ECEF frame to ENU frame
T=          [-sin(r_long)             cos(r_long)            0          0;
             -sin(r_lat)*cos(r_long) -sin(r_lat)*sin(r_long) cos(r_lat) 0;
              cos(r_lat)*cos(r_long)  cos(r_lat)*sin(r_long) sin(r_lat) 0;
              0                       0                      0          1];          

Qp_ENU = T*Qp_ECEF*T';     