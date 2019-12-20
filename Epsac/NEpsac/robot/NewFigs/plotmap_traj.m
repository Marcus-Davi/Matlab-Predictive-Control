% clear;clc
fig = openfig('traj_b2'); 
map = load('map_matlab.mat');
map = map.map_matlab;

MATRIX = map.occupancyMatrix;

axObjs = fig.Children;
dataObjs = axObjs(2).Children

%Robot
x_ref = dataObjs(2).XData;
y_ref = dataObjs(2).YData;


%Ref
x_rob = dataObjs(1).XData;
y_rob = dataObjs(1).YData;
close all
show(map)
hold on
plot(x_rob,y_rob,'red','linewidth',2)
plot(x_ref,y_ref,'--black')
