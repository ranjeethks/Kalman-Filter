%Ranjeeth KS, University of Calgary
function Qr = computeQr(method,Qr_r,Qr_b,number_of_sats);

if(method == 'singleP')
    Qr= Qr_r; 
elseif(method == 'singleD')
    m = eye(number_of_sats);
    B=[m -m];

    temp = [Qr_r zeros(number_of_sats,number_of_sats);
            zeros(number_of_sats,number_of_sats) Qr_b];

    Qr =  B*temp*B';
elseif(method == 'doubleD')
    m1 = ones(number_of_sats-1,1);
    m2 = eye(number_of_sats-1);

    B = [m1 -m2 -m1 m2];

    temp = [Qr_r zeros(number_of_sats,number_of_sats);
            zeros(number_of_sats,number_of_sats) Qr_b];

    Qr =  B*temp*B';  
end