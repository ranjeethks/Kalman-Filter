%Ranjeeth KS, University of Calgary
function [dx] = computeLS(h_constraint,dz, H, R)


dx = (inv(H'*inv(R)*H))*(H'*inv(R)*dz);
