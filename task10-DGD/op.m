function K=op(M,P,B,l,x)
jk=0;
jk2=1;
for k=1:l
    jk2=jk2-M(index(l),index(k));
end
for k=l+1:no_of_tasks
    jl=1;
    for j=1:l-1
        j1=j1-M(index(j),index(k));
    end
    jk=jk+min(x*P(k,l),j1);
end;
jk-jk2;
end