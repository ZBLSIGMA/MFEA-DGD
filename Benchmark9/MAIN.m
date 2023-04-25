clear all
pop_M=100; % population size 100
gen=1000; % generation count 1000
p_il = 0; % probability of individual learning (BFGA quasi-Newton Algorithm) --> Indiviudal Learning is an IMPORTANT component of the MFEA.
rmp=0.7;% random mating probability
reps = 20; % repetitions 20

gama=0.07;
sigma=0.01; %sigma <= 0 indicates randomly selection

xoperator="DGDX";

moperator="DGDM";

%######################
gama=0.07;
sigma=0;
xoperator="SBX";  % SBX
moperator="DGDM";  %PM
for index =1:9
    disp('benchmark')
    disp(index)
    Tasks = benchmark(index);
    %MFEA_GHS_data(index)=MFEA_GHS(Tasks,pop_M,gen,rmp,p_il,reps);  
    %MFEA_GHS_data(index)=MFEAAKT(Tasks,pop_M,gen,rmp,p_il,reps);  
    MFEA_GHS_data(index)=MFEA_DGD(Tasks,pop_M,gen,rmp,p_il,reps,gama, sigma,xoperator,moperator);  
    
end
gama = num2str(gama,1);
sigma = num2str(sigma, 2);
path = strcat('task9_result_DGD(gama=', gama, ',sigma=',sigma,',x=',xoperator,',m=',moperator,').mat');
save(path,'MFEA_GHS_data');


