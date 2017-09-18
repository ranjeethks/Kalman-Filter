%Ranjeeth KS, University of Calgary
function [InpDataStatic InpDataKinematic InpDataSimulator epochSt epochKi epochSim]= ReadData ()

fid1 = fopen('D:\Study\ENGO 620\Lab2\Code\Static.bin' ,'rb');
k = 0; % k stands for k'th epoch
while(k<4000)
    k=k+1; 
    InpDataStatic.gpstime(k) = fread(fid1,1,'double');
    InpDataStatic.numsats(k) = fread(fid1,1,'char');
    for m = 1:InpDataStatic.numsats(k)
        % m stands for m'th satellite
    InpDataStatic.satdata(k).prnID(m) = fread(fid1,1,'char');
    InpDataStatic.satdata(k).sat_data(m,:) = fread(fid1,10,'double');
    end
    
end
epochSt = k;
fclose(fid1);

fid1 = fopen('D:\Study\ENGO 620\Lab2\Code\Kinematic.bin' ,'rb');
k = 0; % k stands for k'th epoch
while(k<1190)
    
    k=k+1; 
    InpDataKinematic.gpstime(k) = fread(fid1,1,'double');
    InpDataKinematic.numsats(k) = fread(fid1,1,'char');
    for m = 1:InpDataKinematic.numsats(k)
        % m stands for m'th satellite
    InpDataKinematic.satdata(k).prnID(m) = fread(fid1,1,'char');
    InpDataKinematic.satdata(k).sat_data(m,:) = fread(fid1,10,'double');
    end
   
    
end
 epochKi = k;
fclose(fid1);
    
fid1 = fopen('D:\Study\ENGO 620\Lab2\Code\Simulator.bin' ,'rb');
k = 0; % k stands for k'th epoch
while(k<650)
    k=k+1; 
    InpDataSimulator.gpstime(k) = fread(fid1,1,'double');
    InpDataSimulator.numsats(k) = fread(fid1,1,'char');
    for m = 1:InpDataSimulator.numsats(k)
        % m stands for m'th satellite
    InpDataSimulator.satdata(k).prnID(m) = fread(fid1,1,'char');
    InpDataSimulator.satdata(k).sat_data(m,:) = fread(fid1,10,'double');
    end
    
end
epochSim = k;
fclose(fid1);
        
