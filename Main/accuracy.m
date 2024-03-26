
pre_label=[ones(50,1);ones(50,1).*2;ones(50,1).*3];
idx = [index1, index2, index3]';
idx = index';
pre_label =label;
[Acc,rand_index,match]=AccMeasure(pre_label,idx);
NMI_val = nmi(pre_label, idx);
NMI_val = NMI_val*100;
