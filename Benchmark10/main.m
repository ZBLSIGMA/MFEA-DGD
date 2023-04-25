clear all

pop_M=100; % population size 100
gen=1000; % generation count 1000
p_il = 0; % probability of individual learning (BFGA quasi-Newton Algorithm) --> Indiviudal Learning is an IMPORTANT component of the MFEA.
rmp=0.7;% random mating probability
reps = 1; % repetitions 20

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gama=0.07;
sigma=0; %sigma <= 0 indicates randomly selection
xoperator="SBX";  % SBX  DGDX
moperator="DGDM";  % PM  DGDM

for index =1:10
%for index =1:10
    Tasks = benchmark(index);
    disp('benchmark')
    disp(index);
    MFEA_DGD_data(index)=MFEA_DGD(Tasks,pop_M,gen,rmp,p_il,reps,gama,sigma,xoperator,moperator);  
end
gama = num2str(gama,1);
sigma = num2str(sigma, 2);
path = strcat('task10_result_DGD(gama=', gama, ',sigma=',sigma,',x=',xoperator,',m=',moperator,').mat');
%'task10_result_DGD.mat'
save(path,'MFEA_DGD_data');
