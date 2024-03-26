function [center, tmpu, dist, obj_fcn] = afkmmean(data, cluster_n, options)

if nargin ~= 2 && nargin ~= 3
error('Too many or too few input arguments!');
end

data_n = size(data, 1);
in_n = size(data, 2);
% Change the following to set default options
default_options = [2;	% exponent for the partition matrix U
100;	% max. number of iteration
1e-5;	% min. amount of improvement
1];	% info display during iteration 

if nargin == 2
options = default_options;
else
% If "options" is not fully specified, pad it with default values.
if length(options) < 4
tmp = default_options;
tmp(1:length(options)) = options;
options = tmp;
end
% If some entries of "options" are nan's, replace them with defaults.
nan_index = find(isnan(options)==1);
options(nan_index) = default_options(nan_index);
if options(1) <= 1
error('The exponent should be greater than 1!');
end
end

expo = options(1);	% Exponent for U
max_iter = options(2);	% Max. iteration
min_impro = options(3);	% Min. improvement
display = options(4);	% Display info or not
lambda =10;  
%% 10 for sead and 500 for iris
obj_fcn = zeros(max_iter, 1);	% Array for objective function

% U = initfkmu(cluster_n, data_n);
% center =0.1*rand(cluster_n, in_n);
b = max(max(data));
a = min(min(data));
center = a + (b-a) * rand(cluster_n, in_n);
% Initial fuzzy partition
% Main loop
for i = 1:max_iter
    %[U, center, dist, obj_fcn(i)] = afkmstep(data, U, cluster_n, center, i);    % with different
   
    outw = distfkmw(center, data);
    dist = distfcm(center, data);  
    
         %% mf update
     for j= 1: cluster_n
           temp = outw{j};
           den(j,:) = exp((-sum(temp,2)/lambda))';
     end

           den_sum = sum(den,1);

      for j= 1: cluster_n
         temp = outw{j};
         for l=1:data_n 
             tmpu(j,l) = exp(-sum(temp(l,:))/lambda)/den_sum(l);
         end 
      end
      
      obj_fcn(i) = sum(sum((dist.^2).*tmpu)) +  sum(lambda* sum(tmpu'.*log(tmpu'),2));
      %====================================================================================
      
      
      center = tmpu*data./(sum(tmpu,2)*ones(1,size(data,2))); %new center
      
       

if display 
fprintf('Iteration count = %d, obj. fcn = %f\n', i, obj_fcn(i));
end
% check termination condition
    if i > 3
        if (abs(obj_fcn(i) == obj_fcn(i-1)))
            return
        end
    end
    
%
%             dist = distfcm(center, data); 
%             obj_fcn(i) = sum(sum((dist.^2).*tmpu)) +  sum(lambda* sum(tmpu'.*log(tmpu'),2)); 
%             if (abs(obj_fcn(i) == obj_fcn(i-1))) 
%                 return
%             end 
%             
%         end
%     end
%     if display 
% fprintf('Iteration count = %d, obj. fcn = %f\n', i, obj_fcn(i));
% end
end

iter_n = i;	% Actual number of iterations 
obj_fcn(iter_n+1:max_iter) = [];