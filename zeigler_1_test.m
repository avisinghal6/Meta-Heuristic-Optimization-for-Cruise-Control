
close all; clear all; clc;

% define pant transfer function
s = tf('s');
%opt = stepDataOptions('StepAmplitude',45);
sys = 2.4767/(s^3+6.0476*s^2+5.2856*s+0.238);
%sys=7.43*exp(-4.95*s)/(1+50.4832*s);
% Obtain step response of the system
[y,t] = step(sys);
plot(t,y,'LineWidth',2); grid on; xlabel('Time(s)'); ylabel('Amplitude');
title('Open loop Response');

%% Obtain Inflection point and Draw Tangent
yp = diff(y);
ypp = diff(y,2);
% Find the root using FZERO
t_infl = fzero(@(T) interp1(t(2:end-1),ypp,T,'linear','extrap'),0);
y_infl = interp1(t,y,t_infl,'linear');
hold on;
plot(t_infl,y_infl,'ro');

% Draw Tangent at inflection point
h = mean(diff(t));
dy = gradient(y, h);
[~,idx] = max(dy);
b = [t([idx-1,idx+1]) ones(2,1)] \ y([idx-1,idx+1]);            % Regression Line Around Maximum Derivative
tv = [-b(2)/b(1); (1-b(2))/b(1)];                               % Independent Variable Range For Tangent Line Plot
f = [tv ones(2,1)] * b;                                         % Calculate Tangent Line

plot(tv, f, '-r','LineWidth',1.5)
ylim([0 max(y)]);

%% finding T and L
L = tv(1);
T = tv(2) - tv(1);

% PID parameters
a = L/T;
Kp = 1.2/a;
Ti = 2*L;
Td = L/2;
Ki=Kp/Ti;
Kd=Kp*Td;

% cont = Kp(1 + 1/(s*Ti) + s*Td);
cont = Kp + Kp/(s*Ti) + Kp*Td*s;

cl_sys = feedback(cont*sys,1);
t = [0:0.01:100];
[yc,tc] = step(cl_sys,t);
figure;
plot(tc,yc,'LineWidth',2); xlabel('Time(s)'); ylabel('Amplitude');
title('Zeigler Nicholas Optimized Closed Loop Response');
grid on;
STI2 = stepinfo(yc,tc,1);
    ST2=STI2.SettlingTime
    PO2=STI2.Overshoot
    RT2 = STI2.RiseTime
    Kp
    Ki
    Kd


