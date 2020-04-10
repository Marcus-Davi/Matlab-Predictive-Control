% clear;clc

fig1 = openfig('traj_b1'); 
fig2 = openfig('traj_b2'); 
map = load('map_matlab.mat');
map = map.map_matlab;

MATRIX = map.occupancyMatrix;

axObjs1 = fig1.Children;
dataObjs1 = axObjs1(2).Children

axObjs2 = fig2.Children;
dataObjs2 = axObjs2(2).Children

% Ref
x_ref = dataObjs1(2).XData;
y_ref = dataObjs1(2).YData;


% EPSAC B1
x_rob1 = dataObjs1(1).XData;
y_rob1 = dataObjs1(1).YData;
% EPSAC B1
x_rob2 = dataObjs2(1).XData;
y_rob2 = dataObjs2(1).YData;

close all
show(map)
hold on
plot(x_ref,y_ref,'--black')
plot(x_rob1,y_rob1,'red','linewidth',2)
plot(x_rob2,y_rob2,'blue:','linewidth',2)
legend('Reference','EPSAC-B_1','EPSAC-B_2')
grid on
