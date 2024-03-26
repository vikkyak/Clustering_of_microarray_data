clear all;
clc;
addpath('./datasets')
load("iris.data")
X = iris(:,1:end-1);
actual_label = iris(:,end);
Nc = 3;  % numer of actual clusters
M = 2;   % fuzzification parameter must be greater than one  
maxIter = 10000;
minImprove = 1e-6;
clusteringOptions = [M maxIter minImprove true];
[gamma, center, U, W, distout] = SCADCul(X,Nc,clusteringOptions);  % this function also have SCAD1 and SCAD2
[maxU , index]= max(U);
[val1, index1] = find(U(1,:) == maxU);
[val2, index2] = find(U(2,:) == maxU);
[val3, index3] = find(U(3,:) == maxU);
param.c=Nc;
param.m=2;
param.e=1e-6;
param.val=1;
data.X = X;
result.data.f=U';
result.data.d=distout';
result.cluster.v=center;
for param = 1:3
result = validity(result,data,param);
end
pridicted_label = index';
[Acc,rand_index,match]=AccMeasure(actual_label,pridicted_label);
