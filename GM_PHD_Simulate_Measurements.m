%GM_PHD_Simulate_Measurements
%Matlab code by Bryan Clarke b.clarke@acfr.usyd.edu.au 

%This file generates simulated measurement data for  the simulation
%described in Vo&Ma.
%There will be gaussian noise on the measurement and Poisson-distributed clutter
%in the environment. 

%If you want to use this PHD filter implementation for another problem, you
%will need to replace this script with another one that populates Z,
%zTrue, and simMeasurementHistory (Z is used in a lot of the update code,
%zTrue and simMeasurementHistory are used in GM_PHD_Simulate_Plot)

%Note: It is possible to get no measurements if the target is not detected
%and there is no clutter
s = sprintf('Step Sim: Simulating measurements.');
disp(s);

%Simulate target movement
simTarget1State = F * simTarget1State;
simTarget2State = F * simTarget2State;

%Save target movement for plotting
simTarget1History = [simTarget1History, simTarget1State];
simTarget2History = [simTarget2History, simTarget2State];


%First, we generate some clutter in the environment.
%clutter = zeros(2,nClutter);%The observations are of the form [x; y]
for i = 1:nClutter
    clutterX =tau(k,i)*100;   
    if i >1 && tau(k,i) == tau(k,i-1)
        break;
    end
    clutter(1,i) = clutterX;
    clutter(2,i) = 0;
end


%Append clutter
Z = clutter;

%Store history
simMeasurementHistory{k} =  Z;