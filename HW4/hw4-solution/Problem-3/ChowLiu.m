clear; 
close all;
import brml.*
load('ChowLiuData.mat');
figure;
drawNet(ChowLiu(X)); title('Chow Liu Net from data')