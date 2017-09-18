%Ranjeeth KS, University of Calgary
function DOPs = computeDOP(Qp_ENU) 
 EDOP=sqrt(Qp_ENU(1,1));
 NDOP=sqrt(Qp_ENU(2,2));
 VDOP=sqrt(Qp_ENU(3,3));
 HDOP=sqrt(EDOP*EDOP+NDOP*NDOP);
 PDOP=sqrt(HDOP*HDOP+VDOP*VDOP);
 
 DOPs=[EDOP NDOP VDOP HDOP PDOP];