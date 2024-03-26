function [delta, tmpu,tmpW, center, dist, obj_fcn] = SCADCulstep(data, U, W, cluster_n, center, expo, delta)
%% cost function

data_n = size(data, 1);
in_n = size(data, 2);
%eta= 100;
q = 1.5;
mf = U.^expo;      % MF matrix after exponential modification
[dist, distw] = distfkm(center, data, W, q);       % fill the distance matrix
outw = distfkmw(center, data);

%% SCAD1
% obj_fcn = sum(sum((dist.^2).*mf)) +  delta'*sum(W*W',2);  
%% SCAD2
% obj_fcn = sum(sum((distw.^2).*mf));  
%% Proposed
 obj_fcn = sum(sum((dist.^2).*mf)) +  delta'*sum((W*(log(W)-ones(cluster_n,in_n))'),2);  

%% weight update for SCAD1

% for j= 1: cluster_n
% 
%     temp = cell2mat(outw(j));
%     for l=1:in_n
%         tmpW(j,l) =  1/in_n  +    ((mf(j,:).^expo)*(((dist(j,:).^2)./in_n)' - temp(:,l)))/(2*delta(j));
%         if tmpW(j,l)<=0
%             tmpW(j,l) = 0;
%         else
%             tmpW(j,l) = tmpW(j,l);
%         end
%     end
% 
% end

%% weight update for SCAD2

% for j= 1: cluster_n   %% for neumerator
%     temp = cell2mat(outw(j));
%     for l=1:in_n
%         num(j,l) = 1/(mf(j,:)*temp(:,l))^(1/(q-1));
%     end
% end
% 
% for j= 1: cluster_n   %% for denomerator
%     temp = cell2mat(outw(j));
%     den(j,:) = 1./(mf(j,:)*temp).^(1/(q-1));
% end
% den_sum = sum(den,2);
% 
% tmpW = num./den_sum;


%% our weight update in clustering to remove trivial solutions in the weight update

  for j = 1:cluster_n
      
      delta(j,:) = 5*abs(((mf(j,:).^expo)*(distw(j,:).^2)')/(W(j,:)*(log(W(j,:))-ones(1, in_n))'));
      
  end 

for j= 1: cluster_n   %% for neumerator
    temp = cell2mat(outw(j));
    for l=1:in_n
        num(j,l) = exp(-((mf(j,:).^expo)*temp(:,l))/delta(j));
    end
end

for j= 1: cluster_n   %% for denomerator
    temp = cell2mat(outw(j));
    den(j,:) = exp((-(mf(j,:).^expo)*temp)/delta(j));
end
den_sum = sum(den,2);

tmpW = num./den_sum;


%% mf update for SCAD1

% for j= 1: cluster_n
%     temp = outw{j};
%     den(j,:) = (1./(tmpW(j,:)*temp')).^(1/(expo -1));
% end
% 
% den_sum = sum(den,1);
% 
% for j= 1: cluster_n
%     temp = outw{j};
%     for i=1:data_n
%         tmpu(j,i) = 1/((tmpW(j,:)*temp(i,:)')^(1/(expo -1))*den_sum(i));
%     end
% end

%% mf update for SCAD2

% for j= 1: cluster_n
%     temp = outw{j};
%     denm(j,:) = (1./((tmpW(j,:).^q)*temp')).^(1/(expo -1));
% end
% 
% den_summ = sum(denm,1);
% 
% for j= 1: cluster_n
%     temp = outw{j};
%     for i=1:data_n
%         tmpu(j,i) = 1/(((tmpW(j,:).^q)*temp(i,:)')^(1/(expo -1))*den_summ(i));
%     end
% end


%% mf update for in proposed

for j= 1: cluster_n
    temp = outw{j};
    denm(j,:) = (1./((tmpW(j,:))*temp')).^(1/(expo -1));
end

den_summ = sum(denm,1);

for j= 1: cluster_n
    temp = outw{j};
    for i=1:data_n
        tmpu(j,i) = 1/(((tmpW(j,:))*temp(i,:)')^(1/(expo -1))*den_summ(i));
    end
end

%% center update

center = ((tmpu.^expo)*data)./(sum(tmpu.^expo,2)*ones(1,size(data,2)));

% delta update for SCAD1

% for j = 1:cluster_n
% 
%     delta(j,:) = 10*(tmpu(j,:)*(distw(j,:))')/(tmpW(j,:)*tmpW(j,:)');
% 
% end

end
