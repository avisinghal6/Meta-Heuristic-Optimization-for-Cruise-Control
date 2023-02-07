
N = 7;                                                                      % number of parameters to optimize 
M = 500;                                                                    % number of particles
T = 5000;                                                                     % iteration number
x_min = [-1 -1 -1 -1 -1.1 -1 -1 ];                                                            % parameter minimum limit
x_max = [1 1 1 1 1 1 1 ];                                                         % parameter maximum limit
delta_min = -0.5;                                                              % Neighborhood minimum limit
delta_max = 0.5;                                                            % Neighborhood maximum limit
t = 1;

%     for n = 1:N,
%         ip(1,n) = x_min(1,n)+(x_max(1,n)-x_min(1,n))*rand(1,1);
%     end
ip=[0.0554,0.1661,0.1661,0.0554,1.0000,-1.0784,0.6450,-0.1236];

theta1=[ip(1),ip(2),ip(3),ip(4)];
theta2=[ip(5),ip(6),ip(7),ip(8)];
[h,w]=freqz(theta1,theta2,2000);
% x_min1=[0 0 0 0 0 0 0 0];
% x_max1=[0 0 0 0 0 0 0 0];
% if(ip(1)>0)
%     x_min1(1)=0;
%     x_max1(1)=1;
% else x_min1(1)=-1;
%     x_max1(1)=0;
%     
% end
% if(ip(2)>0)
%     x_min1(2)=0;
%     x_max1(2)=1;
% else x_min1(2)=-1;
%     x_max1(2)=0;
%     
% end
% if(ip(3)>0)
%     x_min1(3)=0;
%     x_max1(3)=1;
% else x_min1(3)=-1;
%     x_max1(3)=0;
%     
% end
% if(ip(4)>0)
%     x_min1(4)=0;
%     x_max1(4)=1;
% else x_min1(4)=-1;
%     x_max1(4)=0;
%     
% end
% if(ip(5)>0)
%     x_min1(5)=0;
%     x_max1(5)=1;
% else x_min1(5)=-1;
%     x_max1(5)=0;
%     
% end
% if(ip(6)>0)
%     x_min1(6)=0;
%     x_max1(6)=1;
% else x_min1(6)=-1.01;
%     x_max1(6)=0;
%     
% end
% if(ip(7)>0)
%     x_min1(7)=0;
%     x_max1(7)=1;
% else x_min1(7)=-1;
%     x_max1(7)=0;
%     
% end
% if(ip(8)>0)
%     x_min1(8)=0;
%     x_max1(8)=1;
% else x_min1(8)=-1;
%     x_max1(8)=0;
%     
% end

%[h,w]=freqz(theta1,theta2);
%[ph,w]=phasez(theta1,theta2);


for m = 1:M,
    for n = 1:N,
        x(m,n) = x_min(1,n)+(x_max(1,n)-x_min(1,n))*rand(1,1);
    end
end

fx_b = ones(1,M)';  

xf = [x fx_b]; 

for i = 1:M,
    b0 = xf(i,1);
    b1 = xf(i,2);
    b2 = xf(i,3);
    b3 = xf(i,4);
    a0 = 1;
    a1 = xf(i,5);
    a2 = xf(i,6);
    a3 = xf(i,7);
    
    theta3=[b0,b1,b2,b3];
    theta4=[a0,a1,a2,a3];
    [h1,w]=freqz(theta3,theta4,2000);
    %[ph1,w]=phasez(theta3,theta4);
    
    e(i,1)=sum((abs(h)-abs(h1)).^2);
    %e1(i,1)=sum((abs(ph)-abs(ph1)).^2)/512;
    %fe1(i,1)=e(i,1)*(e1(i,1).^2);\

     
    xf(i,8) = e(i,1);
    
    
end

xf_sira = sortrows(xf,[8]);   
x_best = xf_sira(1:M,1:N);                                                  % best order parameter values
f_best_son = xf_sira(1,8); 
cw=zeros(M/2,1);
%f_best=f_best_son;
for lamda =1:M/2
    cw(lamda,1)=xf_sira(lamda,8);
end

g_son = x_best(1,:);                                                        %the solution with the best fit
fbs = 0.2*ones(T-1,1);
while t <= T,                                                               % iteration start
    t=t
    iyi_xf = xf_sira(1:M/2,1:N+1);
            
    for md = 1:M/2,                                                         % plus / minus delta neighborhood
        for nd = 1:N+1,
            if(xf_sira(md,8)>cw)
        delta_min=delta_min-0.01;
        delta_max=delta_max+0.01;
        cw=xf_sira(md,8);
    else if(xf_sira(md,8)<cw)
             delta_min=delta_min+0.01;
             delta_max=delta_max-0.01;
             cw=xf_sira(md,8);
        else if(xf_sira(md,8)==cw)
                delta_min=delta_min;
                delta_max=delta_max;
                cw=xf_sira(md,8);
            end
        end
    end
            delta = delta_min+(delta_max-delta_min)*unifrnd(-1,1);
                komsu_xf(md,nd) = iyi_xf(md,nd)+delta;
           
        end
    end
     for i=1:M/2
         for j=1:N
            if komsu_xf(i,j)<x_min(j)
               komsu_xf(i,j)=x_min(j);
            elseif komsu_xf>x_max(j)
                   komsu_xf=x_max(j);
            end
         end
      end
    %**********************************************************************
    
    for m = 1:M/2,                                                          % the remaining 50% are reproduced
        for n = 1:N,
            x_kalan(m,n) = x_min(1,n)+(x_max(1,n)-x_min(1,n))*rand(1,1);
        end
    end
    yeni_xf = [x_kalan ones(M/2,1)];                                  %the remaining 50% were randomly produced
    xf_sira_yeni = [komsu_xf;yeni_xf];                               % new candidate matrix (including eligibility)
    %**********************************************************************

     for i = 1:M,
        b0 = xf_sira_yeni(i,1);
        b1 = xf_sira_yeni(i,2);
        b2 = xf_sira_yeni(i,3);
        b3 = xf_sira_yeni(i,4);
        a0 = 1;
        a1 = xf_sira_yeni(i,5);
        a2 = xf_sira_yeni(i,6);
        a3 = xf_sira_yeni(i,7);
    
        theta3=[b0,b1,b2,b3];
        theta4=[a0,a1,a2,a3];
       
   
        [h1,w]=freqz(theta3,theta4,2000);
    %[ph1,w]=phasez(theta3,theta4);
        e(i,1)=sum((abs(h)-abs(h1)).^2);
    %e1(i,1)=sum((abs(ph)-abs(ph1)).^2)/512;
    %fe1(i,1)=e(i,1)*(e1(i,1).^2);
    


     
   % xf(i,5) = e(i,1);                                                        % total error
        xf_sira_yeni(i,8) = e(i,1);
    
    
    end
    xf_sira = sortrows(xf_sira_yeni,[8]); 
    
    x_best = xf_sira(1:M,1:N);                                              % en iyi sýralamalý parametre deðerleri
    f_best = xf_sira(1,8);                                                  % en iyi uygunluk
    g = x_best(1,:);                                                        % en iyi uygunluða sahip çözüm
    %**********************************************************************   
    if f_best < f_best_son,                                                 % en iyi uygunluk karþýlaþtýr
        f_best_son = f_best;                                                % en iyi uygunluðu sakla
        g_son = g;                                                          % en iyi sonucu sakla
    end
    t=t+1;                                                                  % iterasyon sayacýný bir arttýr
    f_best_son
    g_son
    %fbs(t) = f_best_son;
end
sonuc = g_son
p1=[g_son(1),g_son(2),g_son(3),g_son(4)];
p2=[1,g_son(5),g_son(6),g_son(7)];

%teta1=subs(theta3,[x,x1],[a,b]);

%double(teta1);











%vpa

        
       
       

