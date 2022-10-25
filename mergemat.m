a1=load('/Users/owner/Desktop/Oct_2022_Imp/new_curlyI_best.mat');
f1=fieldnames(a1)
a2=load('/Users/owner/Desktop/Oct_2022_Imp/new_curlyI.mat');
f2=fieldnames(a2);
v=[a1.(f1{1});a2.(f2{1})]
save /Users/owner/Desktop/Oct_2022_Imp/all_and_best_Oct_2022.mat v