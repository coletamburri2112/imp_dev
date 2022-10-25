%process qpp data from Laura Hayes and compare to impulsiveness index

% 16 September - redo for relative impulsiveness
% 25 October 2022 - redo for new selection of impulsiveness

clear;
qpp_file = '/Users/owner/Desktop/CU_Research/QPP_Study/qpp_hayes_data.csv';
imp_file = '/Users/owner/Desktop/Oct_2022_Imp/imp_dev/all_and_best_Oct_2022.mat';
qpp_data = readtable(qpp_file);
imp_data = load(imp_file,'curly_Is_best','curly_Is_relative_best',...
    'bestflaresname');
    
qpp_start = table2array(qpp_data(:,2));
qpp_peak = table2array(qpp_data(:,3));
qpp_end = table2array(qpp_data(:,4));
qpp_xrcl = table2array(qpp_data(:,5));
qpp_period = table2array(qpp_data(:,6));
qpp_detected = table2array(qpp_data(:,7));
qpp_flareduration = table2array(qpp_data(:,8));
qpp_clcorr = table2array(qpp_data(:,9));

imp_ind = imp_data.curly_Is_relative_best;
imp_name = cell2mat(imp_data.bestflaresname);
bfnames = imp_data.bestflaresname;

qpp_startmat = cellstr(qpp_start);
qpp_peakmat = cellstr(qpp_peak);

bestimp_qpp_date_relative = {};
bestimp_qpp_period_relative = NaN(500,2);
for j=1:500
    impflname = imp_name(j,:);
    imp_dt = [impflname(1:8),impflname(10:13)];
    for i=1:length(qpp_start)
        
        qpp_starti = qpp_startmat(i,:);
        qpp_starti = qpp_starti{1};
        qpp_dt = [qpp_starti(1:4),qpp_starti(6:7),qpp_starti(9:10),qpp_starti(12:13),qpp_starti(15:16)];
        if strcmp(imp_dt,qpp_dt)==1
            bestimp_qpp_date_relative{j}={impflname};
            bestimp_qpp_period_relative(j,1)=qpp_period(i);
            
        end
    end
end

bestimp_qpp_period_relative(:,2) = imp_ind(:,2);
save(imp_file,'bestimp_qpp_period_relative','bestimp_qpp_date_relative','-append');

