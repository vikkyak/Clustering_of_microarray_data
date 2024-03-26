function [center, tmpu, tmpW, dist, obj_fcn] = ewkmmean(data, cluster_n, options)

if nargin ~= 2 && nargin ~= 3
error('Too many or too few input arguments!');
end

data_n = size(data, 1);
in_n = size(data, 2);
fea_m = in_n;
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
lambda =10000;  
%% 10 for sead and 500 for iris
obj_fcn = zeros(max_iter, 1);	% Array for objective function
b = max(max(data));
a = min(min(data));
center = a + (b-a) * rand(cluster_n, in_n);
%center = 1+100*rand(cluster_n, in_n);
tmpW = initfkmw (data,cluster_n,fea_m);
% Initial fuzzy partition
% Main loop

for i = 1:max_iter
    %[U, center, dist, obj_fcn(i)] = afkmstep(data, U, cluster_n, center, i);    % with different
   
    [~, dist] = distewkm(center, data, tmpW);
    
    %% partition matrix update
    [mindst, indmin] = min(dist);
    for j= 1: cluster_n
        for l = 1:data_n
            if (dist(j,l)<= mindst(l))
                tmpu(j,l) = 1;
            else
                tmpu(j,l) = 0;
            end
        end
    end

%  [mindst, indexm] = min(dist);
%  for j= 1: cluster_n
%      for l = 1:data_n
%          for r= 1: cluster_n
%              if (dist(j,l)== dist(r,l))
%                  ind = find(dist(j,l) == dist(r,l));
%                  tmpu(ind,l) = 1;
%              elseif (dist(j,l)<= mindst(l))
%                  tmpu(j,l) = 1;
%              else
%                  tmpu(j,l) = 0;
%              end
%          end
%      end
%  end


% for j= 1: cluster_n
%     for l = 1:data_n
%         for r= 1: cluster_n
%             if dist(j,l)<= dist(r,l)
%                 tmpu(j,l) = 1;
%             else
%                 tmpu(j,l) = 0;
%             end
%         end
%     end
% end



    %% center update
     center = tmpu*data./(sum(tmpu,2)*ones(1,size(data,2)));
     outw = distfkmw(center, data);
     
     %% weight update
     for j= 1: cluster_n
         temp = cell2mat(outw(j));
         for l=1:in_n
             tmpW(j,l) =exp(-(tmpu(j,:)*temp(:,l)/lambda))/sum(exp((-tmpu(j,:)*temp)./lambda));
         end
     end
     %%
    
  obj_fcn(i) = sum(sum((dist.^2).*tmpu)) +  sum(lambda* sum(tmpW'.*log(tmpW'),2));
      %====================================================================================

      if display
          fprintf('Iteration count = %d, obj. fcn = %f\n', i, obj_fcn(i));
      end
% check termination condition
if i > 1
    if abs(obj_fcn(i) - obj_fcn(i-1)) < min_impro, break; end
end
    
end

iter_n = i;	% Actual number of iterations 
obj_fcn(iter_n+1:max_iter) = [];