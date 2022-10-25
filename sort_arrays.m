load('/Users/owner/Desktop/Oct_2022_Imp/imp_dev/all_and_best_Oct_2022.mat')
rsqurs_mult = rsqurs_exp2_four2(:,1).*rsqurs_exp2_four2(:,2);
length(rsqurs_mult);
for i=1:length(rsqurs_mult)
    num(i) = i;
end
r2s_with_num = [rsqurs_mult num'];
sorted = sortrows(r2s_with_num);


for i=1:length(sorted(:,1))
    if ~isnan(sorted(i,1))
        remnansort(i,:) = sorted(i,1:2);
    end
end

flip = flipud(remnansort);
best = flip(1:500,2);
bestsorted = sort(best);

save('/Users/owner/Desktop/Oct_2022_Imp/imp_dev/all_and_best_Oct_2022.mat','bestsorted','-append')