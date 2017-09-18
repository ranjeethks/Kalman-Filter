%Ranjeeth KS, University of Calgary
function [QMatrix] = GetProcessNoiseMatrix(states,dt,spectral_density_vector,number_of_states)

SDV = spectral_density_vector;

if(states == 'PCO')
    
    if(number_of_states == 5)
    QMatrix = [SDV(1)*dt        0           0                  0                         0;
               0            SDV(2)*dt       0                  0                         0;
               0                0       SDV(3)*dt              0                         0;
               0                0           0       SDV(4)*dt + SDV(8)*((dt).^3)/3      SDV(8)*(dt.^2)/2;
               0                0           0           SDV(8)*(dt.^2)/2                SDV(8)*dt];
           
    elseif(number_of_states == 4)
    QMatrix = [SDV(1)*dt        0           0                  0    ;
               0            SDV(2)*dt       0                  0    ;
               0                0       SDV(3)*dt              0    ;
               0                0           0           SDV(4)*dt   ];
      
    end
    
elseif(states == 'PVC')
    
    QMatrix(1,1) = SDV(5)*((dt).^3)/3;
    QMatrix(1,5) = SDV(5)*((dt).^2)/2;
    
    QMatrix(2,2) = SDV(6)*((dt).^3)/3;
    QMatrix(2,6) = SDV(6)*((dt).^2)/2;
    
    QMatrix(3,3) = SDV(7)*((dt).^3)/3;
    QMatrix(3,7) = SDV(7)*((dt).^2)/2;
    
    QMatrix(4,4) = SDV(4)*dt + SDV(8)*((dt).^3)/3;
    QMatrix(4,8) = SDV(8)*((dt).^2)/2;
    
    QMatrix(5,1) = SDV(5)*((dt).^2)/2;
    QMatrix(5,5) = SDV(5)*dt;
    
    QMatrix(6,2) = SDV(6)*((dt).^2)/2;
    QMatrix(6,6) = SDV(6)*dt;
    
    QMatrix(7,3) = SDV(7)*((dt).^2)/2;
    QMatrix(7,7) = SDV(7)*dt;
    
    QMatrix(8,4) = SDV(8)*((dt).^2)/2;
    QMatrix(8,8) = SDV(8)*dt;
    
end


