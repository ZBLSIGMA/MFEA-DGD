classdef Chromosome    
    properties
        rnvec; % (genotype)--> decode to find design variables --> (phenotype) 
        factorial_costs;
        factorial_ranks;
        scalar_fitness;
        skill_factor;
        survival=0;
        p_il;
        options;
        grad;
    end    
    methods        
        function object = initialize(object,D,p_il,options)            
            object.rnvec = rand(1,D);
            object.p_il = p_il;
            object.options = options;
        end
        
        function [object,calls] = evaluate_vec(object,Tasks,p_il,no_of_tasks,options)
            calls = 0;
            if object.skill_factor == 0
                calls=0;
                for i = 1:no_of_tasks
                    [object.factorial_costs(i),xxx,funcCount]=fnceval(Tasks(i),object.rnvec,p_il,options);
                    calls = calls + funcCount;
                end
            else
                object.factorial_costs(1:no_of_tasks)=inf;
                %object.factorial_costs(no_of_tasks)=inf;
                for i = 1:no_of_tasks
                    if object.skill_factor == i
                        [object.factorial_costs(object.skill_factor),object.rnvec,funcCount]=fnceval(Tasks(object.skill_factor),object.rnvec,p_il,options);
                        calls = funcCount;
                        break;
                    end
                end
            end
        end
        
        function [object,calls] = evaluate_SOO(object,Task,p_il,options)   
            [object.factorial_costs,object.rnvec,funcCount]=fnceval(Task,object.rnvec,p_il,options);
            calls = funcCount;
        end
        
        % SBX
        function object=crossover(object,p1,p2,cf)
            object.rnvec=0.5*((1+cf).*p1.rnvec + (1-cf).*p2.rnvec);
            object.rnvec(object.rnvec>1)=1;
            object.rnvec(object.rnvec<0)=0;
        end

        function object=crossover11(object,p1,p2,cf,QWE,L,sigma)             
            p1.rnvec=p1.rnvec-QWE(1,:).*sigma/L;
            p2.rnvec=p2.rnvec-QWE(2,:).*sigma/L;
            object.rnvec=0.5*((1+cf).*p1.rnvec + (1-cf).*p2.rnvec);
            object.rnvec(object.rnvec>1)=1;
            object.rnvec(object.rnvec<0)=0;
        end

        function object=crossover_sbx(object,p1,p2,cf,QWE,L,sigma,isfirst)             
            %p1.rnvec=p1.rnvec-QWE(1,:).*sigma/L;
            %p2.rnvec=p2.rnvec-QWE(2,:).*sigma/L;
            %object.rnvec=0.5*((1+cf).*p1.rnvec + (1-cf).*p2.rnvec);
            %object.rnvec(object.rnvec>1)=1;
            %object.rnvec(object.rnvec<0)=0;
            
            %SBX
            dim=length(p1.rnvec);
            yita1=0.6;
            x_max=1;
            x_min=0;
            u1=zeros(1,dim);
            gama=zeros(1,dim);
            for j=1:dim
               u1(j)=rand(1);
               if u1(j)<0.5
                   gama(j)=(2*u1(j))^(1/(yita1+1));
               else
                   gama(j)=(1/(2*(1-u1(j))))^(1/(yita1+1));
               end
               off_1(j)=0.5*((1+gama(j))*p1.rnvec(j)+(1-gama(j))*p2.rnvec(j));
               off_2(j)=0.5*((1-gama(j))*p1.rnvec(j)+(1+gama(j))*p2.rnvec(j));
               %make offsprings in the defined field
               if(off_1(j)>x_max)
                   off_1(j)=x_max;
               elseif(off_1(j)<x_min)
                   off_1(j)=x_min;
               end
               if(off_2(j)>x_max)
                   off_2(j)=x_max;
               elseif(off_2(j)<x_min)
                   off_2(j)=x_min;
               end
            end
            if isfirst==1
                object.rnvec=off_1;
            else
                object.rnvec=off_2;
            end
            object.rnvec(object.rnvec>1)=1;
            object.rnvec(object.rnvec<0)=0;
        end

        function object=mutate1(object,p,dim,qq,L,sigma)          
            object.rnvec = p.rnvec-qq.*sigma/L;          
        end    

        function object=mutate_pm(object,p,dim,qq,L,sigma)  
            x_num=p.rnvec;
            yita2=0.6;
            x_max=1;
            x_min=0;
            u2=zeros(1,x_num);
            delta=zeros(1,x_num);
            for j=1:x_num
               u2(j)=rand(1);
               if(u2(j)<0.5)
                   delta(j)=(2*u2(j))^(1/(yita2+1))-1;
               else
                   delta(j)=1-(2*(1-u2(j)))^(1/(yita2+1));
               end
               off_1(j)=off_1(j)+delta(j);
               %make offsprings in the defined field
               if(off_1(j)>x_max)
                   off_1(j)=x_max;
               elseif(off_1(j)<x_min)
                   off_1(j)=x_min;
               end
            end
        end
        
        % polynomial mutation
        function object=mutate(object,p,dim,mum)
            rnvec_temp=p.rnvec;
            for i=1:dim
                if rand(1)<1/dim
                    u=rand(1);
                    if u <= 0.5
                        del=(2*u)^(1/(1+mum)) - 1;
                        rnvec_temp(i)=p.rnvec(i) + del*(p.rnvec(i));
                    else
                        del= 1 - (2*(1-u))^(1/(1+mum));
                        rnvec_temp(i)=p.rnvec(i) + del*(1-p.rnvec(i));
                    end
                end
            end  
            object.rnvec = rnvec_temp;          
        end    
    end
end