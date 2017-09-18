%Ranjeeth KS, University of Calgary
function [az el]= calcAzEl(x_ENU)
az = atan2(x_ENU(1),x_ENU(2));
el = asin(x_ENU(3)/sqrt(x_ENU(1).^2+ x_ENU(2).^2 + x_ENU(3).^2));