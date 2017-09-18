%Ranjeeth KS, University of Calgary
function [residue elevation availablePRN MDBout dx_ext_relH dx_ext_relV]=arrangeResidue(check,epochs,MDB1,dx_ext_rel1);
 p=1;
 
 dumPRNarray = 0;
 for k = 1:epochs
     dumPRNarray = [dumPRNarray check(1).PRNIDarray(k).prn];
 end

 %to find what PRNs are present, check upto 35 PRNs
 for prnID=1:35;
     if(find(dumPRNarray==prnID))
         availablePRN(p)=prnID;
         p=p+1;
     end
 end
 
 %residues and elevation angles are written as column elements in 'residue' and 'elevation'
 %PRN number is located in corresponding column of 'availablePRN'
 for k=1:epochs
     if(check(1).num_blunders(k)==0)
         sat_ind=1;
         for n=1:length(availablePRN)
             prnInd = find(check(1).PRNIDarray(k).prn==availablePRN(n));
             if(prnInd)
                 %column-wise sync between 'availablePRN' vector and 'residue' matrix
                 residue(k,sat_ind)=check(1).residual(k).res(prnInd);
                 elevation(k,sat_ind)=check(1).elevationarray(k).ele(prnInd);
             else
                 %IF a PRN is not available for a given epoch, filled with NaN
                 residue(k,sat_ind)=NaN;
                 elevation(k,sat_ind)=NaN;
             end
             sat_ind=sat_ind+1;
         end
     elseif(check(2).num_blunders(k)==0)
         sat_ind=1;
         for n=1:length(availablePRN)
             prnInd = find(check(2).PRNIDarray(k).prn==availablePRN(n));
             if(prnInd)
                 %column-wise sync between 'availablePRN' vector and 'residue' matrix
                 residue(k,sat_ind)=check(2).residual(k).res(prnInd);
                 elevation(k,sat_ind)=check(2).elevationarray(k).ele(prnInd);
             else
                 %IF a PRN is not available for a given epoch, filled with NaN
                 residue(k,sat_ind)=NaN;
                 elevation(k,sat_ind)=NaN;
             end
             sat_ind=sat_ind+1;
         end
     end
 end
 
 for k=1:epochs
     
         sat_ind=1;
         for n=1:length(availablePRN)
             prnInd = find(check(1).PRNIDarray(k).prn==availablePRN(n));
             if(prnInd)
                 %column-wise sync between 'availablePRN' vector and 'residue' matrix
                 MDBout(k,sat_ind)=MDB1(k).sat(prnInd);
                 dx_ext_relH(k,sat_ind)=dx_ext_rel1(k).Nsats(prnInd,1);
                 dx_ext_relV(k,sat_ind)=dx_ext_rel1(k).Nsats(prnInd,2);
             else
                 %IF a PRN is not available for a given epoch, filled with NaN
                 MDBout(k,sat_ind)=NaN;
                 dx_ext_relH(k,sat_ind)=NaN;
                 dx_ext_relV(k,sat_ind)=NaN;
             end
             sat_ind=sat_ind+1;
         end
 end
 
