%Ranjeeth KS, University of Calgary
function [TMatrix] = GetTransitionMatrix(states, dt, number_of_states);

if(states == 'PCO')
    if(number_of_states == 5)
    TMatrix = [1 0 0 0 0;
               0 1 0 0 0;
               0 0 1 0 0;
               0 0 0 1 dt;
               0 0 0 0 1];
    elseif(number_of_states == 4)
    TMatrix = [1 0 0 0 ;
               0 1 0 0 ;
               0 0 1 0 ;
               0 0 0 1 ];
    end

elseif(states == 'PVC')
    TMatrix = [1  0  0  0  dt 0  0  0;
               0  1  0  0  0 dt  0  0;
               0  0  1  0  0  0 dt  0;
               0  0  0  1  0  0  0 dt;
               0  0  0  0  1  0  0  0;
               0  0  0  0  0  1  0  0;
               0  0  0  0  0  0  1  0;
               0  0  0  0  0  0  0  1];
end