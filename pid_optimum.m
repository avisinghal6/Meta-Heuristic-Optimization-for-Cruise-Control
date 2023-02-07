function M= pid_optimum(x)
s=tf('s');
plant= 1.2/(0.00077*s^3+0.0539*s^2+1.441*s);
kp=x(1)
ki=x(2)
kd=x(3)
controller=kp+ki*1/s+kd*s;
%step(feedback(controller*plant,1));
dt=0.01;
t=0:dt:1;
e=1-step(feedback(controller*plant,1),t);
M= sum(t'.*abs(e)*dt);
end
