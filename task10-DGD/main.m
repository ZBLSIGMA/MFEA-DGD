clear all

pop_M=100; % population size 100
gen=1000; % generation count 1000
p_il = 0; % probability of individual learning (BFGA quasi-Newton Algorithm) --> Indiviudal Learning is an IMPORTANT component of the MFEA.
rmp=0.7;% random mating probability
reps = 20; % repetitions 20
parfor index =1:10
%for index =1:10
    Tasks = benchmark(index);
    disp('benchmark')
    disp(index);
    MFEA_DGD_data(index)=MFEA_DGD(Tasks,pop_M,gen,rmp,p_il,reps);  
end
save('task10_result_DGD.mat','MFEA_DGD_data');