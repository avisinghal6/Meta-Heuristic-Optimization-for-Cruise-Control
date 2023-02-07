clear all
close all
clc
%--------------------------------------------------------------------------
p=30; % Population size
c=10; % number of pairs of chromosomes to be crossovered
m=10; % number chromosomes to be mutated
tg=5; % Total number of generations 
%--------------------------------------------------------------------------
figure
title('Blue - Average         Red - Maximum');
xlabel('Generation')
ylabel('Objective Function Value')
hold on
P=population(p);
K=0;
[x1 y1]=size(P);
P1 = 0;
for i=1:tg   
    Cr=crossover(P,c);
    Mu=mutation(P,m);
    P(p+1:p+2*c,:)=Cr;
    P(p+2*c+1:p+2*c+m,:)=Mu;
    E= evaluation(P);
    [P S]=selection(P,E,p);
    K(i,1)=sum(S)/p;
    K(i,2)=S(1);
    plot(K(:,1),'b.'); drawnow
    hold on
    plot(K(:,2),'r.'); drawnow
end
Max_fitness_value=max(K(:,2))
P2=P(1,:);
A=bi2de(P2(1,1:y1/3));
x=3+A*(4-(3))/(2^(y1/3)-1);

B=bi2de(P2(1,y1/3+1:y1));
y=0.1+B*(0.25-(0.1))/((2^y1/3)-1);

C=bi2de(P2(1,y1/3+2:y1));
z=3+B*(4-(3))/((2^y1/3)-1);

Optimal_solution=[x y z]