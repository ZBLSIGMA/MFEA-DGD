clear all

pop_M=100; % population size 100
gen=1000; % generation count 1000
p_il = 0; % probability of individual learning
rmp=0.7;% random mating probability
reps = 20; % repetitions 20

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gama=0.1;
sigma=0;  %sigma <= 0 indicates randomly selection


for index =1:10
    Tasks = benchmark(index);
    disp(['benchmark ', num2str(index)]);
    MFEA_DGD_data(index)=MFEA_DGD(Tasks,pop_M,gen,rmp,p_il,reps,gama,sigma);  
end
gama = num2str(gama,1);
sigma = num2str(sigma, 2);
path = strcat('CEC2021_result_DGD(gama=', gama, ',sigma=',sigma,').mat');
save(path,'MFEA_DGD_data');
