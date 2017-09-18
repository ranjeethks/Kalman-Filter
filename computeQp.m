%Ranjeeth KS, University of Calgary
function [Qp] = computeQp(h_constraint,H,Qr)


Qp = (inv(H'*inv(Qr)*H));

