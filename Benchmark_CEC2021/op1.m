function K=op1(M,E,B,l,x)
jk=0;
jk2=1;
for k=1:l-1
    jk2=jk2-M(index(j),index(l));
end
for j=l+1:no_of_tasks
    jl=1;
    for k=1:l-1
        j1=j1-M(B(j),B(k));
    end
    jk=jk+min(x*E(j,l),j1);
end;
jk-jk2;
end