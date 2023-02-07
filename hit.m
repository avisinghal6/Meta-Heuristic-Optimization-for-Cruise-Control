clc;
clear;
close all

%% Problem Definition
        % Cost Function
nVar=4;             % Number of Dec-ision Variables
VarSize=[1 nVar];   % Decision Variables Matrix Size
VarMin=[3 0 3 0];         % Decision Variables Lower Bound
VarMax=[3 1 4 2500];         % Decision Variables Upper Bound
%% ABC Settings
MaxIt=10;              % Maximum Number of Iterations
nPop=10;               % Population Size (Colony Size)
nOnlooker=nPop;         % Number of Onlooker Bees
L=round(0.6*nVar*nPop); % Abandonment Limit Parameter (Trial Limit)
a=1;                    % Acceleration Coefficient Upper Bound
%% Initialization
% Empty Bee Structure
empty_bee.Position=[];
empty_bee.Cost=[];
% Initialize Population Array
pop=repmat(empty_bee,nPop,1);
% Initialize Best Solution Ever Found
BestSol.Cost=inf;
% Create Initial Population
for i=1:nPop
    for j=1:4
      pop(i).Position(j)=VarMin(1,j)+(VarMax(1,j)-VarMin(1,j))*rand(1,1);
    end
    
     Kp = pop(i).Position(1);
    Ki =pop(i).Position(2);
    Kd = pop(i).Position(3);
   Ni =pop(i).Position(4);
    sim('pidn');                                                     
    load('error1.mat');                                                      
    load('error.mat');                                                      
    e = ITAE(2,:);                                                         
    time = ITAE(1,:);                                                       
    e1=fark(2,:);
    time1=fark(1,:);
    STI = stepinfo(e1,time1,0);
    ST=STI.SettlingTime;
    PO=STI.Peak;
    SSE = e(end);
    pop(i).Cost=SSE;
    if pop(i).Cost<=BestSol.Cost
        BestSol=pop(i);
    end
end
% Abandonment Counter
C=zeros(nPop,1);
% Array to Hold Best Cost Values
BestCost=zeros(MaxIt,1);
%% ABC Main Loop
for it=1:MaxIt
    
    % Recruited Bees
    for i=1:nPop
        
        % Choose k randomly, not equal to i
        K=[1:i-1 i+1:nPop];
        k=K(randi([1 numel(K)]));
        
        % Define Acceleration Coeff.
        phi=a*unifrnd(-1,1,VarSize);
        %phi=a*unifrnd(-1,1);
        %lamda=rand();
        % New Bee Position
        newbee.Position=pop(i).Position+phi.*(pop(i).Position-pop(k).Position);
        
            for b=1:4
               if newbee.Position(b)<VarMin(b)
                   newbee.Position(b)=VarMin(b);
               elseif newbee.Position(b)>VarMax(b)
                  newbee.Position(b)=VarMax(b);
               end
            end
         
        % Evaluation
       
     Kp = newbee.Position(1);
    Ki =newbee.Position(2);
    Kd = newbee.Position(3);
   Ni = newbee.Position(4);
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
    SSE = e(end);
        newbee.Cost= SSE;
        
        % Comparision
        if newbee.Cost<=pop(i).Cost
            pop(i)=newbee;
        else
            C(i)=C(i)+1;
        end
        
    end
    
    % Calculate Fitness Values and Selection Probabilities
    F=zeros(nPop,1);
    MeanCost = mean([pop.Cost]);
    for i=1:nPop
        F(i) = exp(-pop(i).Cost/MeanCost); % Convert Cost to Fitness
    end
    P=F/sum(F);
    
    % Onlooker Bees
    for m=1:nOnlooker
        
        % Select Source Site
        i=RouletteWheelSelection(P);
        
        % Choose k randomly, not equal to i
        K=[1:i-1 i+1:nPop];
        k=K(randi([1 numel(K)]));
        
        % Define Acceleration Coeff.
        phi=a*unifrnd(-1,1,VarSize);
        
        % New Bee Position
        newbee.Position=pop(i).Position+phi.*(pop(i).Position-pop(k).Position);
        for b=1:4
               if newbee.Position(b)<VarMin(b)
                   newbee.Position(b)=VarMin(b);
               elseif newbee.Position(b)>VarMax(b)
                  newbee.Position(b)=VarMax(b);
               end
            end
        
        % Evaluation
         Kp = newbee.Position(1);
    Ki =newbee.Position(2);
    Kd = newbee.Position(3);
    Ni = newbee.Position(4);
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
    SSE = e(end);
        newbee.Cost= SSE;
        
        % Comparision
        if newbee.Cost<=pop(i).Cost
            pop(i)=newbee;
        else
            C(i)=C(i)+1;
        end
        
    end
    
    % Scout Bees
    for i=1:nPop
        if C(i)>=L
            for j=1:4
             pop(i).Position(j)=VarMin(1,j)+(VarMax(1,j)-VarMin(1,j))*rand(1,1);
            end
    Kp = pop(i).Position(1);
    Ki =pop(i).Position(2);
    Kd = pop(i).Position(3);
    Ni = pop(i).Position(4);
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
    SSE = e(end);
            pop(i).Cost=SSE;
            C(i)=0;
        end
    end
    
    % Update Best Solution Ever Found
    for i=1:nPop
        if pop(i).Cost<=BestSol.Cost
            BestSol=pop(i);
        end
    end
    
    % Store Best Cost Ever Found
    BestCost(it)=BestSol.Cost;
    
    % Display Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    
end
    
%% Results
figure;
%plot(BestCost,'LineWidth',2);
semilogy(BestCost,'LineWidth',2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;
Kp=BestSol.Position(1)
Ki=BestSol.Position(2)
Kd=BestSol.Position(3)
Ni=BestSol.Position(4);
sim('pidn');

  load('error2.mat');

    e2 = fark2(2,:);
    time2 = fark2(1,:);
    STI2 = stepinfo(e2,time2,1);
    ST2=STI2.SettlingTime
    PO2=STI2.Overshoot
    RT2 = STI2.RiseTime
