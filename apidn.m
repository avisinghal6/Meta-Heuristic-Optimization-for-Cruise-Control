clear all
close all
clc

%-------------------------------------------------------------------------
p=25; % Population size
c=10; % number of pairs of chromosomes to be crossovered
m=10; % number chromosomes to be mutated
tg=10; % Total number of generations 
%--------------------------------------------------------------------------
figure
title('Blue - Average         Red - Maximum');
xlabel('Generation')
ylabel('Objective Function Value')
hold on
P=round(rand(p,40));
K=0;
[x1,y1]=size(P);
P1 = 0;
for o=1:tg
    Z=zeros(2*c,y1); 
for i = 1:c
    r1=randi(x1,1,2);
    while r1(1)==r1(2)
        r1=randi(x1,1,2);
    end
    A1=P(r1(1),:); % parent 1
    A2=P(r1(2),:); % parent 2
    r2=1+randi(y1-1);
    B1=A1(1,r2:y1);
    A1(1,r2:y1)=A2(1,r2:40);
    A2(1,r2:40)=B1;
    Z(2*i-1,:)=A1;
    Z(2*i,:)=A2;
end
Cr=Z;
W=zeros(m,y1);
for i = 1:m
    r1=randi(x1);
    A1=P(r1,:); % random parent
    r2=randi(y1);
    if A1(1,r2)== 1
        A1(1,r2) = 0; % flick the bit
    else
        A1(1,r2) = 1;
    end
    W(i,:)=A1;
end
Mu=W;
 P(p+1:p+2*c,:)=Cr;
    P(p+2*c+1:p+2*c+m,:)=Mu;
    H=zeros(1,x1);
for i = 1:x1
    A=bi2de(P(i,1:y1/4));
    Kp=18+A*(20-(18))/(2^(y1/4)-1);

    B=bi2de(P(i,y1/4+1:2*(y1/4)));
    Ki=B*(1-(0))/(2^(y1/3)-1);

    C=bi2de(P(i,2*(y1/4)+1:3*(y1/4)));
    Kd=C*(1-(0))/(2^(y1/4)-1);
    
     NEW=bi2de(P(i,3*(y1/4)+1:y1));
    Ni=NEW*(800-(0))/(2^(y1/4)-1);
 sim('pidn');                                                     % proses simülasyonu
    load('error1.mat');                                                      % prosesten hata deðerleri çekiliyor
    load('error.mat');                                                      % prosesten hata deðerleri çekiliyor
    e = ITAE(2,:);                                                          % error vector
    time = ITAE(1,:);                                                       % time vector
    e1=fark(2,:);
    time1=fark(1,:);
    STI = stepinfo(e1,time1,0);
    ST=STI.SettlingTime;
    PO=STI.Peak;
    %SSE = (0.75*ST)+(0.18*e(end))+(0.07*PO); 
    SSE=e(end);
    
    H(1,i)= 1/(SSE+1);
end
Y1 = zeros(p,y1);
%F = F + 10; % adding 10 to ensure no chromosome has negative fitness
% elite selection
elite=3;
for i = 1:elite
    [r1,c1]=find(H==max(H));
    Y1(i,:)=P(max(c1),:); 
    P(max(c1),:)=[];
    Fn(i)=H(max(c1));
    H(:,max(c1))=[];
end
D=H/sum(H);
E=cumsum(D);
N=rand(1);
d1=1;
d2=elite;
while d2 <=p-elite
    if N<=E(d1)
        Y1(d2+1,:)=P(d1,:);
        Fn(d2+1)=H(d1);
        N=rand(1);
        d2=d2+1;
        d1=1;
    else
        d1=d1+1;
    end
end
P=Y1;
S=Fn;
K(o,1)=sum(S)/p;
    K(o,2)=S(1);
    plot(K(:,1),'b.'); drawnow
    hold on
    plot(K(:,2),'r.'); drawnow
end
Max_fitness_value=max(K(:,2))
P2=P(1,:);
A=bi2de(P2(1,1:y1/4));
    x=18+A*(20-(18))/(2^(y1/4)-1);

    B=bi2de(P2(1,y1/4+1:2*(y1/4)));
    y=B*(1-(0))/(2^(y1/3)-1);

    C=bi2de(P2(1,2*(y1/4)+1:3*(y1/4)));
    z=C*(1-(0))/(2^(y1/4)-1);
    
     NEW=bi2de(P2(1,3*(y1/4)+1:y1));
    x_z=NEW*(800-(0))/(2^(y1/4)-1);

Optimal_solution=[x y z x_z]
Kp=x;
Ki=y;
Kd=z;
Ni=x_z;
sim('pidn');

  load('error2.mat');

    e2 = fark2(2,:);
    time2 = fark2(1,:);
    STI2 = stepinfo(e2,time2,1);
    ST2=STI2.SettlingTime
    PO2=STI2.Overshoot
    RT2 = STI2.RiseTime
    
    
    
    
    