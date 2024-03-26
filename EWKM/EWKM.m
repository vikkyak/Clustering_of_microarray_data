% L. Jing, M. K. Ng, and J. Z. Huang, "An entropy weighting K-means algorithm for subspace clustering of high-dimensional sparsedata,” IEEE Trans. Knowl. Data Eng., vol. 19, no. 8, pp. 1026–1041,Aug. 2007.

clear all;
clc;
addpath('./datasets')
load("iris.data")
X = iris(:,1:end-1);
actual_label = iris(:,end);
Nc = 3;  % numer of actual clusters
M = 2;   % fuzzification parameter must be greater than one  
maxIter = 100;
minImprove = 1e-6;
clusteringOptions = [M maxIter minImprove true];
[center, U, tmpW, distout] = ewkmmean(X,Nc,clusteringOptions);  % agglomerative fuzzy k menas
[maxU , index]= max(U);
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