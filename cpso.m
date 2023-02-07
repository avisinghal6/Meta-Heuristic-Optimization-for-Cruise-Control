%tic
clc
clear all
close all

%ip=1:500;
%rng default
LB=[3 0.1 3 0]; %lower bounds of variables
UB=[4 0.25 4 200]; %upper bounds of variables
% pso parameters values
m=4; % number of variables
n=20; % population size
wmax=0.8; % inertia weight
wmin=0.7; % inertia weight
c1=2; % acceleration factor
c2=2; % acceleration factor
% pso main program----------------------------------------------------start
maxite=10; % set maximum number of iteration 
maxrun=1; % set maximum number of runs need to be
%Ni=100;

for run=1:maxrun
 run
 % pso initialization----------------------------------------------start
 for i=1:n
 for j=1:m
 x0(i,j)=LB(j)+rand(1,1)*(UB(j)-LB(j));
 end
 end
 x=x0; % initial population
 v=0.1*x0; 
 for i=1:n
      Kp = x0(i,1);
        Ki = x0(i,2);
        Kd = x0(i,3);
      Ni = x0(i,4);
        sim('cpidn');                                                 % proses simülasyonu
        load('error1.mat');                                                  % prosesten hata deðerleri çekiliyor                                          
         load('error.mat');
        e = ITAE(2,:);                                                          % error vector
    time = ITAE(1,:);                                                       % time vector
    e1=fark(2,:);
    time1=fark(1,:);
    STI = stepinfo(e1,time1,0);
    ST=STI.SettlingTime;
    PO=STI.Peak;
    SSE = e(end);

 f0(i,1)=SSE;
 end
 [fmin0,index0]=min(f0);
 pbest=x0; % initial pbest
 gbest=x0(index0,:); % initial gbest
 % pso initialization------------------------------------------------end

 % pso algorithm---------------------------------------------------start
 ite=1;
 tolerance=1;
 while ite<=maxite 
%&& tolerance>10^-12
 w=wmax-(wmax-wmin)*ite/maxite; % update inertial weight
 % pso velocity updates
     for i=1:n
     for j=1:m
 v(i,j)=w*v(i,j)+c1*rand()*(pbest(i,j)-x(i,j))...
  +c2*rand()*(gbest(1,j)-x(i,j));
     end
    end
 % pso position update
 for i=1:n
 for j=1:m
 x(i,j)=x(i,j)+v(i,j);
 end
 end
 % handling boundary violations
 for i=1:n
 for j=1:m
 if x(i,j)<LB(j)
 x(i,j)=LB(j)+rand();
 elseif x(i,j)>UB(j)
 x(i,j)=UB(j)-rand();
 end
 end
 end
 % evaluating fitness
 for i=1:n
   Kp = x(i,1);
        Ki = x(i,2);
        Kd = x(i,3);
      Ni = x(i,4);
        sim('cpidn');                                                 % proses simülasyonu
        load('error1.mat');                                                  % prosesten hata deðerleri çekiliyor                                          
         load('error.mat');
        e = ITAE(2,:);                                                          % error vector
    time = ITAE(1,:);                                                       % time vector
    e1=fark(2,:);
    time1=fark(1,:);
    STI = stepinfo(e1,time1,0);
    ST=STI.SettlingTime;
    PO=STI.Peak;
    SSE =e(end);  

 f(i,1)=SSE;
 end
 % updating pbest and fitness
 for i=1:n
 if f(i,1)<f0(i,1)
     pbest(i,:)=x(i,:);
 f0(i,1)=f(i,1);
 end
 end
 [fmin,index]=min(f0); % finding out the best particle
 ffmin(ite,run)=fmin; % storing best fitness
 ffite(run)=ite; % storing iteration count
 % updating gbest and best fitness
 if fmin<fmin0
 gbest=pbest(index,:);
 fmin0=fmin;
 end
 % calculating tolerance
%  if ite>100;
%  tolerance=abs(ffmin(ite-100,run)-fmin0);
%  end
 % displaying iterative results
 if ite==1
 disp(sprintf('Iteration Best particle Objective fun'));
 end
 disp(sprintf('%8g %8g %8.4f',ite,index,fmin0));
 ite=ite+1;
 end
 % pso algorithm-----------------------------------------------------end
 gbest;
 Kp = gbest(1);
    Ki = gbest(2);
    Kd = gbest(3);
 Ni=gbest(4);

sim('cpidn');                                                     % proses simülasyonu
    load('error1.mat');                                                      % prosesten hata deðerleri çekiliyor
    load('error.mat');                                                      % prosesten hata deðerleri çekiliyor
    e = ITAE(2,:);                                                          % error vector
    time = ITAE(1,:);                                                       % time vector
    e1=fark(2,:);
    time1=fark(1,:);
    STI = stepinfo(e1,time1,0);
    ST=STI.SettlingTime;
    PO=STI.Peak;
    SSE = e(end);
 fvalue=SSE;
 fff(run)=fvalue;
 rgbest(run,:)=gbest;
 disp(sprintf('--------------------------------------'));
end
% pso main program------------------------------------------------------end
disp(sprintf('\n'));
disp(sprintf('*********************************************************'));
disp(sprintf('Final Results-----------------------------'));
[bestfun,bestrun]=min(fff)
best_variables=rgbest(bestrun,:)
disp(sprintf('*********************************************************'));
Kp = best_variables(1);
    Ki = best_variables(2);
    Kd = best_variables(3);
  Ni=best_variables(4);

sim('cpidn');                                                     % proses simülasyonu
    load('error1.mat');
  load('error.mat');
  load('error2.mat');
e = ITAE(2,:);                                                          % error vector
    time = ITAE(1,:);                                                       % time vector
    e1=fark(2,:);
    time1=fark(1,:);
    STI = stepinfo(e1,time1,0);
    ST=STI.SettlingTime;
    PO=STI.Overshoot;
    e2 = fark2(2,:);
    time2 = fark2(1,:);
    STI2 = stepinfo(e2,time2,1);
    ST2=STI2.SettlingTime
    PO2=STI2.Overshoot
    RT2 = STI2.RiseTime

%toc
% PSO convergence characteristic
% plot(ffmin(1:ffite(bestrun),bestrun),'-k');
% xlabel('Iteration');
% ylabel('Fitness function value');
% title('PSO convergence characteristic')
% %##########################################################################
