%script to make comparison between impulsiveness and other quantities
clear;
imp_file = '/Users/coletamburri/Desktop/imp_dev/all_and_best_Sep_2023.mat';
imp_data = load(imp_file,'curly_Is_best','curly_Is','curly_Is_relative',...
    'curly_Is_relative_best', ...
    'starttimes_corr','maxtimes_corr','endtimes_corr', ...
    'bestflaresname');

rec_data = restore_idl('/Users/coletamburri/Desktop/imp_dev/recratesidl.sav');

flarenames = imp_data.bestflaresname;
pat1 = ["C"]; 
pat2 = ["M"];
pat3 = ["X"];
pat = ["C","M"];
pat_in = contains(flarenames,pat);
pat_in1 = contains(flarenames,pat1);
pat_in2 = contains(flarenames,pat2);
pat_in3 = contains(flarenames,pat3);


imp_ind = imp_data.curly_Is_relative_best;
imp_start = imp_data.starttimes_corr;
imp_peak = imp_data.maxtimes_corr;
imp_end = imp_data.endtimes_corr;
% qpp_per = imp_data.bestimp_qpp_period_relative;
imp_ind_all = imp_data.curly_Is_relative;




recmed = rec_data.MEDIANRECS;
recmean = rec_data.MEANRECS;
recstd = rec_data.STDRECS;
recmaximp = rec_data.MAXIMPRECS;
recmax = rec_data.MAXRECRATES;

% posup = recstd(:,1);
% posdown = recstd(:,1);
% negup = recstd(:,2);
% negdown = recstd(:,2);

% qpp_peron = qpp_per(:,1);
curly_I = imp_ind(:,1);
rise_dur = (imp_peak(:,1)-imp_start(:,1))*24*3600;
dec_dur = (imp_end(:,1)-imp_peak(:,1))*24*3600;
overall_dur = (imp_end(:,1)-imp_start(:,1))*24*3600;
curly_I_all = imp_ind_all(:,1);

curly_I_low = curly_I(pat_in);
% recmean_low = recmean(pat_in,:);
% recmed_low = recmed(pat_in,:);
% recmaximp_low = recmaximp(pat_in,:);
% recmax_low = recmaximp(pat_in,:);

curly_I_low1 = curly_I(pat_in1,:);
% recmean_low1 = recmean(pat_in1,:);
% recmed_low1 = recmed(pat_in1,:);
% recmaximp_low1 = recmaximp(pat_in1,:);
% recmax_low1 = recmaximp(pat_in1,:);

curly_I_low2 = curly_I(pat_in2,:);
% recmean_low2 = recmean(pat_in2,:);
% recmed_low2 = recmed(pat_in2,:);
% recmaximp_low2 = recmaximp(pat_in2,:);
% recmax_low2 = recmaximp(pat_in2,:);

curly_I_low3 = curly_I(pat_in3,:);
% recmean_low3 = recmean(pat_in3,:);
% recmed_low3 = recmed(pat_in3,:);
% recmaximp_low3 = recmaximp(pat_in3,:);
% recmax_low3 = recmaximp(pat_in3,:);

%max rec rates vs. max rec rates in impulsive phase?
% m = figure(7)
% clf
% set(gcf,'Position',[300 300 700 700])
% scatter(recmaximp(:,1),recmax(:,1),25,'fill')
% grid on
% pbaspect([1 1 1])
% xlabel('Max Pos. Rec. Rate, Impulsive Phase [Mx/s]','Fontsize',15);
% ylabel('Max Pos. Rec. rate, Overall [Mx/s]','Fontsize',15);
% title('Maximum Reconnection Rate Comparisons','Fontsize',20);
% saveas(m,'/Users/owner/Desktop/rec_comp_4May/maxreccomp.fig');
% saveas(m,'/Users/owner/Desktop/rec_comp_4May/maxreccomp.png');
% 
%comparison to max rec rate overall
%all flares
l=figure(5)
clf
colormap cool
%subplot(2,1,1)
set(gcf,'Position',[300 300 550 550])
idx = isfinite(recmax(:,1)) & isfinite(curly_I);
fitted0 = fit(recmax(idx,1)/1e18,curly_I(idx),'m*x+b'); 

k=plot(fitted0,'red');
hold on
set(k,'LineWidth',2)
scatter(recmax(:,1)/1e18,curly_I,200,curly_I,'.');
grid on
xlabel('\textsf{Reconnection Rate} [$$10^{18} Mx/s$$]','fontsize',20,'interpreter','latex')
ylabel('\textsf{Impulsiveness Index} [$$ln(min^{-1})$$]','fontsize',20,'interpreter','latex')
ax=gca;
ax.FontSize = 20;
xlim([0 35])
legend off
a=colorbar;
hl = ylabel(a,'\textsf{Impulsiveness}','fontsize',20,'interpreter','latex');

% 
title('$$i$$ \textsf{vs. Reconnection Rate,} $$r^2 = 0.029$$','fontsize',25,'interpreter','latex')

% 
% % pbaspect([1 1 1])
% % 
% % subplot(2,1,2)
% set(gcf,'Position',[300 300 1300 500])
% scatter(curly_I,recmaximp(:,2)/1e18,60,'k.');
% ylabel('Neg. Max. Rec. Rate During Impulsive Phase [10^{18} Mx/s]','fontsize',12)
% xlabel('Impulsiveness [1/s]','fontsize',12)
% title('Impulsiveness vs. Max. Rec. Rate in Impulsive Phase','fontsize',15)
% % % 
% pbaspect([1 1 1])
% 
% saveas(l,'/Users/owner/Desktop/Research/MAT_SOURCE/impcomp/imp_maximprec_comp_allclass.fig');
% saveas(l,'/Users/owner/Desktop/Research/MAT_SOURCE/impcomp/imp_maximprec_comp_allclass.png');
% saveas(l,'/Users/owner/Desktop/fig9.png');
%C/M class
% l=figure(5)
% clf
% subplot(1,2,1)
% set(gcf,'Position',[300 300 1300 500])
% scatter(curly_I_low,recmaximp_low(:,1),60,'k.');
% 
% ylabel('Pos. Max. Rec. Rate During Impulsive Phase [Mx/s]','fontsize',12)
% xlabel('Impulsiveness [mW/m^2*s]','fontsize',12)
% title('Impulsiveness vs. Max. Rec. Rate in Impulsive Phase, C/M Class','fontsize',15)
% 
% pbaspect([1 1 1])
% 
% subplot(1,2,2)
% set(gcf,'Position',[300 300 1300 500])
% scatter(curly_I_low,recmaximp_low(:,2),60,'k.');
% ylabel('Neg. Max. Rec. Rate During Impulsive Phase [Mx/s]','fontsize',12)
% xlabel('Impulsiveness [mW/m^2*s]','fontsize',12)
% title('Impulsiveness vs. Max. Rec. Rate in Impulsive Phase, C/M Class','fontsize',15)
% 
% pbaspect([1 1 1])
% 
% saveas(l,'/Users/owner/Desktop/imp_maximprec_comp_c_and_m.fig');
% saveas(l,'/Users/owner/Desktop/imp_maximprec_comp_c_and_m.png');
% 
% %C class
% l=figure(5)
% clf
% subplot(1,2,1)
% set(gcf,'Position',[300 300 1300 500])
% scatter(curly_I_low1,recmaximp_low1(:,1),60,'k.');
% 
% ylabel('Pos. Max. Rec. Rate During Impulsive Phase [Mx/s]','fontsize',12)
% xlabel('Impulsiveness [mW/m^2*s]','fontsize',12)
% title('Impulsiveness vs. Max. Rec. Rate in Impulsive Phase, C Class','fontsize',15)
% 
% pbaspect([1 1 1])
% 
% subplot(1,2,2)
% set(gcf,'Position',[300 300 1300 500])
% scatter(curly_I_low1,recmaximp_low1(:,2),60,'k.');
% ylabel('Neg. Max. Rec. Rate During Impulsive Phase [Mx/s]','fontsize',12)
% xlabel('Impulsiveness [mW/m^2*s]','fontsize',12)
% title('Impulsiveness vs. Max. Rec. Rate in Impulsive Phase, C Class','fontsize',15)
% 
% pbaspect([1 1 1])
% 
% saveas(l,'/Users/owner/Desktop/imp_maximprec_comp_c.fig');
% saveas(l,'/Users/owner/Desktop/imp_maximprec_comp_c.png');
% 
% %M class
% l=figure(5)
% clf
% subplot(1,2,1)
% set(gcf,'Position',[300 300 1300 500])
% scatter(curly_I_low2,recmaximp_low2(:,1),60,'k.');
% 
% ylabel('Pos. Max. Rec. Rate During Impulsive Phase [Mx/s]','fontsize',12)
% xlabel('Impulsiveness [mW/m^2*s]','fontsize',12)
% title('Impulsiveness vs. Max. Rec. Rate in Impulsive Phase, M Class','fontsize',15)
% 
% pbaspect([1 1 1])
% 
% subplot(1,2,2)
% set(gcf,'Position',[300 300 1300 500])
% scatter(curly_I_low2,recmaximp_low2(:,2),60,'k.');
% ylabel('Neg. Max. Rec. Rate During Impulsive Phase [Mx/s]','fontsize',12)
% xlabel('Impulsiveness [mW/m^2*s]','fontsize',12)
% title('Impulsiveness vs. Max. Rec. Rate in Impulsive Phase, M Class','fontsize',15)
% 
% pbaspect([1 1 1])
% 
% saveas(l,'/Users/owner/Desktop/imp_maximprec_comp_m.fig');
% saveas(l,'/Users/owner/Desktop/imp_maximprec_comp_m.png');
% 
% %X class
% l=figure(5)
% clf
% subplot(1,2,1)
% set(gcf,'Position',[300 300 1300 500])
% scatter(curly_I_low3,recmaximp_low3(:,1),60,'k.');
% 
% ylabel('Pos. Max. Rec. Rate During Impulsive Phase [Mx/s]','fontsize',12)
% xlabel('Impulsiveness [mW/m^2*s]','fontsize',12)
% title('Impulsiveness vs. Max. Rec. Rate in Impulsive Phase, X Class','fontsize',15)
% 
% pbaspect([1 1 1])
% 
% subplot(1,2,2)
% set(gcf,'Position',[300 300 1300 500])
% scatter(curly_I_low3,recmaximp_low3(:,2),60,'k.');
% ylabel('Neg. Max. Rec. Rate During Impulsive Phase [Mx/s]','fontsize',12)
% xlabel('Impulsiveness [mW/m^2*s]','fontsize',12)
% title('Impulsiveness vs. Max. Rec. Rate in Impulsive Phase, X Class','fontsize',15)
% 
% pbaspect([1 1 1])
% 
% saveas(l,'/Users/owner/Desktop/imp_maximprec_comp_x.fig');
% saveas(l,'/Users/owner/Desktop/imp_maximprec_comp_x.png');
% 

%identify bad flare periods
t=1;
for i = 1:length(curly_I)
    
    if isnan(curly_I(i))
        t=t;
    else
        rise_dur_nonan(t)=rise_dur(i);
        dec_dur_nonan(t)=dec_dur(i);  
        t=t+1;
    end
end



%REC RATE AND I
% l=figure(5)
% clf
% subplot(1,2,1)
% set(gcf,'Position',[300 300 1000 500])
% scatter(curly_I_low,recmed_low(:,1),60,'k.');
% 
% ylabel('Pos. Median Reconnection Rate [Mx/s]','fontsize',12)
% xlabel('Impulsiveness [mW/m^2*s]','fontsize',12)
% title('Impulsiveness vs. Pos. Median Rec. Rate, C/M Class','fontsize',15)
% 
% pbaspect([1 1 1])
% 
% subplot(1,2,2)
% set(gcf,'Position',[300 300 1000 500])
% scatter(curly_I_low,recmed_low(:,2),60,'k.');
% ylabel('Neg. Median Reconnection Rate [Mx/s]','fontsize',12)
% xlabel('Impulsiveness [mW/m^2*s]','fontsize',12)
% title('Impulsiveness vs. Neg. Median Rec. Rate, C/M Class','fontsize',15)
% 
% pbaspect([1 1 1])
% 
% saveas(l,'/Users/owner/Desktop/Research/MAT_SOURCE/impcomp/imp_medianrec_comp_c_and_m.fig');
% saveas(l,'/Users/owner/Desktop/Research/MAT_SOURCE/impcomp/imp_medianrec_comp_c_and_m.png');
% 
% q=figure(6)
% clf
% subplot(1,2,1)
% set(gcf,'Position',[300 300 1000 500])
% scatter(curly_I_low,recmean_low(:,1),60,'k.');
% ylabel('Pos. Mean Reconnection Rate [Mx/s]','fontsize',12)
% xlabel('Impulsiveness [mW/m^2*s]','fontsize',12)
% title('Impulsiveness vs. Pos. Mean Rec. Rate, C/M Class','fontsize',15)
% 
% pbaspect([1 1 1])
% 
% subplot(1,2,2)
% set(gcf,'Position',[300 300 1000 500])
% scatter(curly_I_low,recmean_low(:,2),60,'k.');
% ylabel('Neg. Mean Reconnection Rate [Mx/s]','fontsize',12)
% xlabel('Impulsiveness [mW/m^2*s]','fontsize',12)
% title('Impulsiveness vs. Neg. Mean Rec. Rate, C/M class','fontsize',15)
% 
% pbaspect([1 1 1])
% 
% saveas(q,'/Users/owner/Desktop/Research/MAT_SOURCE/impcomp/imp_meanrec_comp_c_and_m.fig');
% saveas(q,'/Users/owner/Desktop/Research/MAT_SOURCE/impcomp/imp_meanrec_comp_c_and_m.png');
% 
% n=figure(7)
% clf
% subplot(1,2,1)
% set(gcf,'Position',[300 300 1000 500])
% scatter(recmed_low(:,2),recmed_low(:,1),60,'k.');
% ylabel('Pos. Median Reconnection Rate [Mx/s]','fontsize',12)
% xlabel('Neg. Median Reconnection Rate [Mx/s]','fontsize',12)
% title('Pos. and Neg. Median Rec. Rate, C/M Class','fontsize',15)
% 
% pbaspect([1 1 1])
% 
% subplot(1,2,2)
% set(gcf,'Position',[300 300 1000 500])
% scatter(recmean_low(:,2),recmean_low(:,1),60,'k.');
% ylabel('Pos. Mean Reconnection Rate [Mx/s]','fontsize',12)
% xlabel('Neg. Mean Reconnection Rate [Mx/s]','fontsize',12)
% title('Pos. and Neg. Mean Reconnection Rate, C/M Class','fontsize',15)
% 
% pbaspect([1 1 1])
% 
% saveas(n,'/Users/owner/Desktop/Research/MAT_SOURCE/impcomp/posnegrecrate_c_and_m.fig');
% saveas(n,'/Users/owner/Desktop/Research/MAT_SOURCE/impcomp/posnegrecrate_c_and_m.png');

% % % % QPP PERIOD AND I
% n=figure(1);
% clf
% colormap cool
% set(gcf,'Position',[300 300 500 500])
% idx2 = (log(qpp_peron) > 4.5);
% qpp_peron(idx2)
% qpp_peron(idx2) = NaN;
% scatter(log(qpp_peron/60),curly_I,200,curly_I,'.');
% idx = isfinite(curly_I) & isfinite(qpp_peron);
% fitted2 = fit(log(qpp_peron(idx)/60),curly_I(idx),'m*x+b'); 
% 
% 
% 
% 
% 
% title('$$i$$ \textsf{vs. QPP Period} ($$r^2 = 0.225$$)','FontSize',28,'interpreter','latex');
% hold on
% p=plot(fitted2,'red');
% set(p,'LineWidth',2)
% ylabel('\textsf{Impulsiveness Index} [$$ln(min^{-1})$$]','fontsize',20,'interpreter','latex');
% xlabel('\textsf{QPP Period} [$$ln(min)$$]','fontsize',20,'interpreter','latex');
% dim = [.6 .55 .3 .3];
% str = {['m = ',num2str(fitted2.m)],['b = ',num2str(fitted2.b)],['r^2 = 0.225']};
% ax=gca;
% ax.FontSize=20;
% %annotation('textbox',dim,'String',str,'FitBoxToText','on','fontsize',14);
% grid on
% legend off
% pbaspect([1 1 1])
% a=colorbar;
% hl = ylabel(a,'\textsf{Impulsiveness}','fontsize',20,'interpreter','latex');

% saveas(n,'/Users/owner/Desktop/Research/MAT_SOURCE/impcomp/imp_qppper_more.fig');
% saveas(n,'/Users/owner/Desktop/Research/MAT_SOURCE/impcomp/imp_qppper_more.png');
% saveas(n,'/Users/owner/Desktop/qppcomp.png');


% incidence of QPP period with impulsiveness, using all flares for bins
% high res bins
% numbin=150;
% 
% 
% edges = min(curly_I_all):((max(curly_I_all))-min(curly_I_all))/numbin:max(curly_I_all);
% 
% 
% notnan = ~isnan(qpp_peron);
% inan = isnan(qpp_peron);
% imp_qpp_notnan=curly_I;
% 
% for i=1:500
%     if inan(i)==1
%         0
%         imp_qpp_notnan(i)=NaN
%     end
% end
% 
% imp_noqpp = curly_I(inan);
% 
% m=figure(2);
% clf
% set(gca,'yscale','log')
% set(gcf,'Position',[300 300 800 500]);
% g=histogram(curly_I,'binedges',edges,'facecolor','magenta','facealpha',1);
% %xlim([-2.5 2.5])
% hold on
% 
% h=histogram(imp_qpp_notnan,numbin,'binedges',edges,'facecolor','blue','facealpha',1);
% legend('All Events','Events with QPPs','fontsize',20,'interpreter','latex')
% grid on;
% xlabel('Impulsiveness [$$ln(min^{-1})$$]','fontsize',20,'interpreter','latex');
% ylabel('Count','fontsize',25,'interpreter','latex');
% 
% %set(gca,'fontsize',15)
% title('Incidence of QPPs in Impulsiveness Distribution, High Res Bins','fontsize',25,'interpreter','latex');
% 
% % saveas(f,'/Users/owner/Desktop/imp_qpp_hist_more.fig');
% % saveas(f,'/Users/owner/Desktop/imp_qpp_hist_more.png');
% 
% all=hist(curly_I,edges);
% part=hist(imp_qpp_notnan,edges);
% noqpppart=hist(imp_noqpp,edges);
% frac = part./all;
% noqppfrac = noqpppart./all;
% fracrel = part/sum(part);
% noqpprel = noqpppart/sum(noqpppart);
% 
% %low res bins
% numbin=14;
% 
% 
% edges2 = min(curly_I_all):((max(curly_I_all))-min(curly_I_all))/numbin:max(curly_I_all);
% 
% 
% notnan = ~isnan(qpp_peron);
% inan = isnan(qpp_peron);
% imp_qpp_notnan=curly_I;
% 
% for i=1:500
%     if inan(i)==1
%         0
%         imp_qpp_notnan(i)=NaN
%     end
% end
% 
% imp_noqpp = curly_I(inan);
% 
% m=figure(3);
% clf
% set(gca,'yscale','log')
% set(gcf,'Position',[300 300 800 500]);
% g2=histogram(curly_I,'binedges',edges2,'facecolor','magenta','facealpha',1);
% %xlim([-2.5 2.5])
% hold on
% 
% h2=histogram(imp_qpp_notnan,numbin,'binedges',edges2,'facecolor','blue','facealpha',1);
% legend('All Events','Events with QPPs','fontsize',20,'interpreter','latex')
% grid on;
% xlabel('Impulsiveness [$$ln(min^{-1})$$]','fontsize',20,'interpreter','latex');
% ylabel('Count','fontsize',25,'interpreter','latex');
% 
% %set(gca,'fontsize',15)
% title('Incidence of QPPs in Impulsiveness Distribution','fontsize',25,'interpreter','latex');
% xlim([-9 0])
% % saveas(f,'/Users/owner/Desktop/imp_qpp_hist_more.fig');
% % saveas(f,'/Users/owner/Desktop/imp_qpp_hist_more.png');
% 
% all2=hist(curly_I,edges2);
% part2=hist(imp_qpp_notnan,edges2);
% noqpppart2=hist(imp_noqpp,edges2);
% frac2 = part2./all2;
% noqppfrac2 = noqpppart2./all2;
% fracrel2 = part2/sum(part2);
% noqpprel2 = noqpppart2/sum(noqpppart2);

%conduct two-sample K-S test to compare the sample with QPPs to the sample
%without QPPs
% l=figure(5);
% set(gcf,'Position',[300 300 800 500]);
% m=histogram(imp_noqpp,numbin,'facecolor','black');
% grid on;
% xlabel('Impulsiveness [ln(nW/m^2/s)]','fontsize',14);
% ylabel('Counts','fontsize',14);
% title('Incidence of QPPs in Impulsiveness Bins','fontsize',25);
% 
% % assume error bars will just be the standard error of the estimate, when
% % compared to the non-qpp 500 flare counts
% 
% errors = (abs(fracrel-noqpprel))./sqrt(part);
% 
% errlow = frac-errors;
% errhigh = frac+errors;

% o=figure(4);
% set(gcf,'Position',[300 300 800 500]);
% bar(edges,frac,'facecolor','red')
% grid on
% xlabel('Impulsiveness [ln(nW/m^2/s)]','fontsize',14);
% ylabel('Fraction','fontsize',14);
% title('Fraction of Events with QPPs','fontsize',25);
% 
% hold on
% 
% er = errorbar(edges,frac,errors,errors);    
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  
% 
% hold off

% saveas(o,'/Users/owner/Desktop/frac_imp_qpp_more.fig');
% saveas(o,'/Users/owner/Desktop/frac_imp_qpp_more.png');
% 
% 
% 
% x1 = h.Values;
% x2 = m.Values;
% 
% [decision,pval,ks2stat]=kstest2(x1,x2);
% 
% 
% saveas(l,'/Users/owner/Desktop/imp_noqpp_hist_more.fig');
% saveas(l,'/Users/owner/Desktop/imp_noqpp_hist_more.png');

% parula=fake_parula();
% magma=magma();
% inferno=inferno();
% plasma=plasma();
% viridis=viridis();
% %%comparison of rise to decay phase
% l=figure(5)
% set(gcf,'Position',[300 300 1300 500])
% 
% clf
% subplot(1,12,1:5)
% rise_dur_nonan(403)=[];
% dec_dur_nonan(403)=[]; 
% curly_I_cp=curly_I;
% 
% rise_dur_cp = rise_dur;
% dec_dur_cp = dec_dur;
% overall_dur_cp = overall_dur;
% curly_I_cp(403)=[];
% rise_dur_cp(403)=[];
% dec_dur_cp(403)=[];
% overall_dur_cp(403)=[];
% idx = isfinite(curly_I_cp) & isfinite(rise_dur_cp);
% %l=figure(5)
% %set(gcf,'Position',[300 300 600 600])
% 
% colormap viridis
% 
% scatter(dec_dur_nonan/60,rise_dur_nonan/60,200,curly_I_cp(idx),'.');
% ax = gca;
% ax.FontSize=15;
% ylabel('Rise Phase Duration $[min]$','fontsize',20,'interpreter','latex')
% xlabel('Decay Phase Duration $[min]$','fontsize',20,'interpreter','latex')
% title('$t_{rise}$ vs. $t_{decay}$','fontsize',30,'interpreter','latex')
% 
% 
% p = polyfit(dec_dur_nonan/60,rise_dur_nonan/60,1);
% yfit = polyval(p,dec_dur_nonan/60);
% hold on
% grid on
% k = plot(dec_dur_nonan/60,yfit,'red');
% set(k,'LineWidth',2)
% 
% a=colorbar;
% hl = ylabel(a,'Impulsiveness','fontsize',20,'interpreter','latex');
% annotation('textbox', [0.31, 0.75, 0.08, 0.05], 'String', '$r^2 = 0.126$','interpreter','latex',fontsize=20)
% 
% subplot(1,12,7:12)
% 
% 
% 
% scatter(rise_dur_cp/60,curly_I_cp,20,'MarkerEdgeColor','#CC6677','MarkerFaceColor','#CC6677');
% hold on
% scatter(dec_dur_cp/60,curly_I_cp,20,'MarkerEdgeColor','#DDCC77','MarkerFaceColor','#DDCC77');
% scatter(overall_dur_cp/60,curly_I_cp,20,'MarkerEdgeColor','#117733','MarkerFaceColor','#117733');
% 
% grid on
% legend('$$t_{rise}$$ ($$r^2 = 0.036$$)','$$t_{dec}$$ ($$r^2 = 0.116$$)','$$t_{flare}$$ ($$r^2 = 0.146$$)','interpreter','latex','Fontsize',20)
% ax = gca;
% ax.FontSize=15;
% idx = isfinite(curly_I_cp) & isfinite(rise_dur_cp);
% 
% fitted0 = fit(log(rise_dur_cp(idx)/60),curly_I_cp(idx),'m*x+b'); 
% 
% fitted1 = fit(log(dec_dur_cp(idx)/60),curly_I_cp(idx),'m*x+b'); 
% fitted2 = fit(log(overall_dur_cp(idx)/60),curly_I_cp(idx),'m*x+b');
% crise = corrcoef(log(rise_dur_cp(idx)/60),curly_I_cp(idx));
% cdec = corrcoef(log(dec_dur_cp(idx)/60),curly_I_cp(idx));
% cover = corrcoef(log(overall_dur_cp(idx)/60),curly_I_cp(idx));
% %k=plot(fitted0,'red');
% %set(k,'LineWidth',2)
% 
% set(gca,'xscale','log')
% ylabel('Impulsiveness','fontsize',20,'interpreter','latex');
% xlabel('Duration [min]','fontsize',20,'interpreter','latex');
% title(['$$i$$ vs. $$t_{rise}$$, $$t_{decay}$$, $$t_{flare}$$'],'FontSize',30,'interpreter','latex');
% xlim([5,500])
% saveas(l,'/Users/owner/Desktop/Research/MAT_SOURCE/impcomp/rise_dec_comp_more.fig');
% saveas(l,'/Users/owner/Desktop/Research/MAT_SOURCE/impcomp/rise_dec_comp_more.png');
% % % 

% DURATIONS AND I
% g=figure(2);
% clf
% set(gcf,'Position',[300 300 1200 500])
% subplot(1,3,1)
% 
% 
% 
% scatter(log(rise_dur/60),curly_I,200,curly_I,'.');
% grid on
% hold on
% 
% idx = isfinite(curly_I) & isfinite(rise_dur);
% pbaspect([1 1 1])
% fitted0 = fit(log(rise_dur(idx)/60),curly_I(idx),'m*x+b'); 
% %k=plot(fitted0,'red');
% %set(k,'LineWidth',2)
% legend off
% ylabel('\textsf{Impulsiveness Index} [$$ln(min^{-1})$$]','fontsize',14,'interpreter','latex');
% xlabel('\textsf{Rise Phase Duration} [$$ln(min)$$]','fontsize',14,'interpreter','latex');
% 
% 
% % dim = [.7 .55 .6 .6];
% % % str = {['m = ',num2str(fitted0.m)],['b = ',num2str(fitted0.b)]};
% % % annotation('textbox',dim,'String',str,'FitBoxToText','on');
% % 
% % saveas(g,'/Users/owner/Desktop/imp_risedur_more_loglog.fig');
% % saveas(g,'/Users/owner/Desktop/imp_risedur_more_loglog.png');
% % saveas(g,'/Users/owner/Desktop/fig6.png');
% % % 
% % h=figure(3);
% % clf
% colormap cool
% ax=gca;
% ax.FontSize=15;
% title(['$$i$$ \textsf{vs.} $$t_{rise}$$ ($$r^2 = 0.053$$)'],'FontSize',22,'interpreter','latex');
% 
% subplot(1,3,2)
% set(gcf,'Position',[300 300 1200 500])
% scatter(log(dec_dur/60),curly_I,200,curly_I,'.');
% grid on
% hold on
% pbaspect([1 1 1])
% 
% %dim = [.7 .55 .6 .6];
% 
% idx = isfinite(curly_I) & isfinite(dec_dur);
% 
% fitted1 = fit(log(dec_dur(idx)/60),curly_I(idx),'m*x+b'); 
% %o=plot(fitted1,'red');
% %set(o,'LineWidth',2)
% str = {['m = ',num2str(fitted1.m)],['b = ',num2str(fitted1.b)]};
% %annotation('textbox',dim,'String',str,'FitBoxToText','on');
% legend off
% ylabel('\textsf{Impulsiveness Index} [$$ln(min^{-1})$$]','fontsize',14,'interpreter','latex');
% xlabel('\textsf{Decay Phase Duration} [$$ln(min)$$]','fontsize',14,'interpreter','latex');

% itle('$$i$$ vs. QPP Period ($$r^2 = 0.288$$)','FontSize',22,'interpreter','latex');
% 
% saveas(h,'/Users/owner/Desktop/imp_decaydur_more_loglog.fig');
% saveas(h,'/Users/owner/Desktop/imp_decaydur_more_loglog.png');
% saveas(h,'/Users/owner/Desktop/fig7.png');
% % 
% l=figure(4);
% clf
% ax=gca;
% ax.FontSize=15;
% title('$$i$$ \textsf{vs.} $$t_{decay}$$ ($$r^2 = 0.083$$)','FontSize',22,'interpreter','latex');
% subplot(1,3,3)
% set(gcf,'Position',[300 300 1200 500])
% scatter(log(overall_dur/60),curly_I,200,curly_I,'.');
% 
% grid on
% hold on
% pbaspect([1 1 1])
% idx = isfinite(curly_I) & isfinite(overall_dur);
% 
% fitted2 = fit(log(overall_dur(idx)/60),curly_I(idx),'m*x+b'); 
% p=plot(fitted2,'red');
% set(p,'LineWidth',2)
% 
% %dim = [.7 .55 .6 .6];
% str = {['m = ',num2str(fitted2.m)],['b = ',num2str(fitted2.b)],['r^2 = 0.1052']};
% %annotation('textbox',dim,'String',str,'FitBoxToText','on','fontsize',12);
% legend off
% ylabel('\textsf{Impulsiveness Index} [$$ln(min^{-1})$$]','fontsize',14,'interpreter','latex');
% xlabel('\textsf{Flare Duration} [$$ln(min)$$]','fontsize',14,'interpreter','latex');
% 
% ax=gca;
% ax.FontSize=15;
% hp4 = get(subplot(1,3,3),'Position')
% a=colorbar('Position', [hp4(1)+.23  hp4(2)+.155  0.01 hp4(4)*0.625])
% hl = ylabel(a,'\textsf{Impulsiveness}','fontsize',20,'interpreter','latex');
% title(['$$i$$ \textsf{vs.} $$t_{flare}$$ ($$r^2 = 0.132$$)'],'FontSize',22,'interpreter','latex');
% % % 
% % 
% % saveas(l,'/Users/owner/Desktop/imp_overdur_more_loglog.fig');
% % saveas(l,'/Users/owner/Desktop/imp_overdur_more_loglog.png');
% % 
% % saveas(l,'/Users/owner/Desktop/fig5.png');