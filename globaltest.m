%Ranjeeth KS, University of Calgary
%%%%%global test%%%%%
function  global_test = globaltest(r,store_R,X2_1_minus_alpha_by_2,X2_alpha_by_2)
 %set test statistics
 zeta = r'*inv(store_R)*r;
 
 if(X2_alpha_by_2> zeta & zeta > X2_1_minus_alpha_by_2)
     global_test = 'pass';
 else
     global_test = 'fail';
 end