% Oscillation of damped pendulum
function [thetadot] = odefun_pendulum(t, theta0, L,m,b);
g = 9.81; % acceleration due to gravity in m/sec^2 
thetadot = zeros(2,1);
thetadot(1) = theta0(2);
thetadot(2) = (-g/L)*sin(theta0(1)) - (b/m)* theta0(2);
end

% Solving 2nd ODE that describes the motion of the pendulum with damping
clear all;close all;clc;
 
t =0:0.1:20; % given time for which pendulum swing in seconds
theta0 = [0 3]; % given intial conditions:@t=0,theta=0rad, dtheta/dt =3rad/sec 
L = 1; % given length of the pendulum in meters
m = 1; % given mass of the pendulum in kg
b = 0.05; % given damping coefficient
 
% Solution for the given ODE
[t, theta] = ode45(@(t,theta) odefun_pendulum(t,theta,L,m,b),t,theta0);
 
% Plotting Angular displacement vs Time
subplot(3,1,1) 
plot(t, theta(:,1))
xlabel('Time in sec')
ylabel('Angular displacement in rad')
 
% Plotting Angular velocity vs Time
subplot(3,1,2)
plot(t, theta(:,2))
xlabel('Time in sec')
ylabel('Angular velocity in rad/sec')
 
% Plotting Angular displacement vs Angular velocity 
subplot(3,1,3)
plot(theta(:,2), theta(:,1))
xlabel('Angular velocity in rad/sec')
ylabel('Angular displacement in rad')
 
% Pendulum Animation
figure(2)
plot(0,0, 'color','r','marker','.') % making origin as the fixed point
xlim([-1 1])  % x-limit of the plot area
ylim([-1.5 0.5])  % y-limit of the plot area
hold on
 
% Plotting pendulum anglar displacement w.r.t space
for i = 1:length(theta) % for all angular displacements
    x = L*sin(theta(i,1)); % x-position of the pendulum
    y = L*-cos(theta(i,1)); % y-position of the pendulum
 P = line([0,x],[0,y]); % creates pendulum position
 bob = viscircles([x y],0.01,'linewidth',25); % creates bob for each pendulum postion
 Q(i) = getframe(gcf); % creates an array of frames of pendulum motion
 pause(0.00001) % execute each frame at specified time
 
 % to delete the previous frame after each execution
 if i < length(theta) 
     delete(P);
     delete(bob);
 end
end
% To create animation
 movie(Q)
 R = VideoWriter('Oscillation of Pendulum','MPEG-4');
 open(R)
 writeVideo(R,Q)
 close(R)
