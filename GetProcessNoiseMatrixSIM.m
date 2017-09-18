%Ranjeeth KS, University of Calgary
function [QMatrix] = GetProcessNoiseMatrix(states,dt,spectral_density_vector,number_of_states, th, K)

SDV = spectral_density_vector;

%following is for 7 states only (xe, xn, bias, ve, vn, drift, thetaDot)
    

    
    QMatrix(1,1) = SDV(4)*((dt).^3)/3 + SDV(7)*K*K*cos(th)*cos(th)*((dt).^5)/20;
    QMatrix(1,2) = SDV(7)*K*K*cos(th)*sin(th)*((dt).^5)/20;
    QMatrix(1,3) = 0;
    QMatrix(1,4) = SDV(4)*((dt).^2)/2 + SDV(7)*K*K*cos(th)*cos(th)*((dt).^4)/8;
    QMatrix(1,5) = SDV(7)*K*K*cos(th)*sin(th)*((dt).^4)/8;
    QMatrix(1,6) = 0;
    QMatrix(1,7) = -SDV(7)*K*cos(th)*((dt).^3)/6;
      
    QMatrix(2,1) = SDV(7)*K*K*cos(th)*sin(th)*((dt).^5)/20;
    QMatrix(2,2) = SDV(5)*((dt).^3)/3 + SDV(7)*K*K*sin(th)*sin(th)*((dt).^5)/20;
    QMatrix(2,3) = 0;
    QMatrix(2,4) = SDV(7)*K*K*sin(th)*cos(th)*((dt).^3)/6;
    QMatrix(2,5) = SDV(5)*((dt).^2)/2 + SDV(7)*K*K*sin(th)*sin(th)*((dt).^4)/8;
    QMatrix(2,6) = 0;
    QMatrix(2,7) = -SDV(7)*K*sin(th)*((dt).^3)/6;
    
    QMatrix(3,1) = 0;
    QMatrix(3,2) = 0;
    QMatrix(3,3) = SDV(3)*(dt) + SDV(6)*((dt).^3)/3;
    QMatrix(3,4) = 0;
    QMatrix(3,5) = 0;
    QMatrix(3,6) = SDV(6)*((dt).^2)/2;
    QMatrix(3,7) = 0;
    
    QMatrix(4,1) = SDV(4)*((dt).^2)/2 + SDV(7)*K*K*cos(th)*cos(th)*((dt).^4)/8;
    QMatrix(4,2) = SDV(7)*K*K*cos(th)*sin(th)*((dt).^4)/8;
    QMatrix(4,3) = 0;
    QMatrix(4,4) = SDV(4)*(dt) + SDV(7)*K*K*cos(th)*cos(th)*((dt).^3)/6;
    QMatrix(4,5) = SDV(7)*K*K*cos(th)*sin(th)*((dt).^3)/6;
    QMatrix(4,6) = 0;
    QMatrix(4,7) = -SDV(7)*K*cos(th)*((dt).^2)/2;
    
    QMatrix(5,1) = SDV(7)*K*K*cos(th)*sin(th)*((dt).^4)/8;
    QMatrix(5,2) = SDV(5)*((dt).^2)/2 + SDV(7)*K*K*sin(th)*sin(th)*((dt).^4)/8;
    QMatrix(5,3) = 0;
    QMatrix(5,4) = SDV(7)*K*K*sin(th)*cos(th)*((dt).^3)/3;
    QMatrix(5,5) = SDV(5)*(dt) + SDV(7)*K*K*sin(th)*sin(th)*((dt).^3)/3;
    QMatrix(5,6) = 0;
    QMatrix(5,7) = -SDV(7)*K*sin(th)*((dt).^2)/2;
    
    QMatrix(6,1) = 0;
    QMatrix(6,2) = 0;
    QMatrix(6,3) = SDV(6)*((dt).^2)/2;
    QMatrix(6,4) = 0;
    QMatrix(6,5) = 0;
    QMatrix(6,6) = SDV(6)*(dt);
    QMatrix(6,7) = 0;
    
    QMatrix(7,1) = -SDV(7)*K*cos(th)*((dt).^3)/6;
    QMatrix(7,2) = -SDV(7)*K*sin(th)*((dt).^3)/6;
    QMatrix(7,3) = 0;
    QMatrix(7,4) = -SDV(7)*K*cos(th)*((dt).^2)/2;
    QMatrix(7,5) = -SDV(7)*K*sin(th)*((dt).^2)/2;
    QMatrix(7,6) = 0;
    QMatrix(7,7) = SDV(7)*(dt);
