%Ranjeeth KS, University of Calgary
function [TMatrix] = GetTransitionMatrixSIM(states, dt, number_of_states, th, K)

%following is for 7 states only (xe, xn, bias, ve, vn, drift, thetaDot)
    TMatrix = [1    0   0   dt      0   0   -K*cos(th)*((dt).^2)/2  ;
               0    1   0   0       dt  0   -K*sin(th)*((dt).^2)/2  ;
               0    0   1   0       0   dt  0                       ;
               0    0   0   1       0   0   -K*cos(th)*(dt)         ;
               0    0   0   0       1   0   -K*sin(th)*(dt)         ;
               0    0   0   0       0   1   0                       ;
               0    0   0   0       0   0   1                      ];

% %following is for 7 states only (xe, xn, bias, ve, vn, drift, thetaDot)
%     TMatrix = [1    0   0   dt      0   0   0  ;
%                0    1   0   0       dt  0   0  ;
%                0    0   1   0       0   dt  0                       ;
%                0    0   0   1       0   0   -K*cos(th)*(dt)         ;
%                0    0   0   0       1   0   -K*sin(th)*(dt)         ;
%                0    0   0   0       0   1   0                       ;
%                0    0   0   0       0   0   1                      ];
