%Cole Tamburri, 2020, University of Colorado at Boulder APS Department
%
%Performs impulsiveness calculations for solar flares in the 304A line as
%recorded by SDO/EVE MEGS-B Level 2 data.  The input data from this
%instrument are downloaded using "downloadsdoevefit.m" as .fit files, and 
%then processed for use by calculating solar quiet values with fft304.m
%(very simple, only one line, which filters the data over 5 hours in order
%to remove flare signatures).  Flares are identified according to the
%ribbondb database (http://solarmuri.ssl.berkeley.edu/~kazachenko/RibbonDB/)
%and the script can be changed to study any of the 2049 flares within this
%database.  

%Start time, peak irradiance, and end time values for the flares are
%calculated in order to obtain values for peak irradiance and full width at
%half height in time.  In addition, models are fit to the rise and decay
%periods of each flare in order to improve the temporal resolution of the
%light curve and obtain a more functional value for FWHH in time (without
%this, the particularly high-impulsiveness flares would be subject to error
%as a result of the large differences in irradiance value between
%consecutive time steps, leading to error in FWHH

%An effort is also made to filter out and correct for low SNR flares.

%clear workspace
 clear;
filename=['/Users/coletamburri/Desktop/Old Mac/CU_Research/sdoeve304_2.mat'];
ribbondb_info = readtable(['/Users/coletamburri/Desktop/imp_dev/ribbondb_v1.0.csv']);
filename2='/Users/coletamburri/Desktop/imp_dev/all_and_best_Sep_2023.mat';
load(filename)
load(filename2)

mkdir('/Users/coletamburri/Desktop/bestflares/')
is = bestsorted;
l=1;

set(gcf,'Position',[100 100 1500 1500])
endj=NaN;
tst=0;
starttimes=datenum(starttimes);
endtimes=datenum(endtimes);
maxtimes =datenum(peaktimes);
risebest =NaN(500,100000);
decaybest =NaN(500,100000);
risetime =NaN(500,100000);
decaytime =NaN(500,100000);
risemod =NaN(500,100000);
decaymod =NaN(500,100000);
curly_Is_best=NaN(500,2);
curly_Is_relative_best=NaN(500,2);
bestflaresname =cell(500,1);
event_curves =NaN(500,100000);
event_times =NaN(500,100000);
starttimes_corr =NaN(500,2);
maxtimes_corr =NaN(500,2);
endtimes_corr =NaN(500,2);

%after verification with ribbondb database, write array which actually
%contains the true indices of the correct flares in the database
vettedbest_corrected =NaN(500,1);
% 
t=0;
fg=1;

for i=1:2049
    event = ribbondb_info{i,2};
    eventstr = event{1};
    in=i;
    

    if i==1957 || i==1967 || i==1977
        event = ribbondb_info{i-1,2};
        eventstr = event{1};
        in=i-1;
    end
    
    if i==735 
        event = ribbondb_info{i+1,2};
        eventstr = event{1};
        in=i+1;
    end
      
         if ismember(i,is)
                t=t+1;
                vettedbest_corrected (t)=in;
                bestflaresname (t) = event;
    

        sqerrstd=[];
        sqerrstd2=[];
        sqerrstd3=[];
        %"windowstart2" is the beginning of the window as defined by ribbondb
        %define our window to be one hour before "windowstart2" to 5 hours after
        %that point (essentially just expand window, having downloaded all
        %SDO/EVE data
        wst=windowstart2(i);
        wstn=windowstart2(i+1);
        %wst - (1.5/24 works for flares 1-25)
        wstb=wst-(2/24);
        wend=wst+(3.5/24);
        
        %find the indices of times in vectime (from SDO/EVE) corresponding to
        %this window
        eventi=find(vectime>wstb & vectime<wend);
        eventi_halfsolarrot=find(vectime>wst-(12/24) & vectime<wend+(12/24));
        
        %assume flares do not overlap
        if wend>wstn
            eventi=find(vectime>wstb & vectime<wstn);
            eventi_halfsolarrot=find(vectime>wst-(12/24) & vectime<(wstn+(12/24)));
        end
        %extract irradiance/time/sq/difference between irradiance and sq
        %values for this event in our defined window
        
        irrev=vec304(eventi);
        
        inan=find(isnan(irrev));
            timeev=vectime(eventi);
    
        sqev=filt304(eventi);
        errev=vecerr(eventi);
        
        passage = vec304(eventi_halfsolarrot);
        
        passaget = vectime(eventi_halfsolarrot);
        %standard deviations of both solar quiet and irradiance values,
        %ignoring NaN
        sqstd=nanstd(sqev);
        irrstd=nanstd(irrev);
        %irrstd=nanstd(passage);
        irrmean=nanmean(passage);
        diff=irrev-sqev;
     
        for p=2:length(irrev)-1
            if irrev(p)>irrev(p-1)+2*irrstd && irrev(p)>irrev(p+1)+2*irrstd 
                irrev(p)=NaN;
            end
        end
     
        for p=2:length(irrev)-1
            if irrev(p)<irrev(p-1)-irrstd  
                irrev(p)=NaN;
            end
        end    
    
              
                
        %N.B.: use solar quiet standard deviation to avoid picking up later flares
        %which might be included in the same window (use irrsq, and later,
        %larger flares might skew the std so that the first, less strong flare
        %might be ignored by the calculation)
        

        
        nans=isnan(irrev);
        
            
        %now begin logic to find the start time of the flare
        %(1) ideal start time scenario - find 40 points which lie consecutively
        %above the solar quiet by 6 times the solar quiet standard deviation;
        %if this is sufficiently far into the window, subtract 80 from that
        %index to find the start time. If this is not sufficiently far into the
        %window, simply subtract 39 and the start time will be roughly the
        %beginning of the window.
        %(2) and (3) perform the same process, but lower the criteria to 30
        %points consecutively above 6*sqstd in (2) and 20 points consecutively
        %above 6*sqstd.  Always subtract an extra 40 points when moving back to
        %find the calculated start time
        %
        n=0;
        m=0;
        lastj=0;
        tst=0;
        for j=1:length(irrev)
            if diff(j)>(2*irrstd)
                n=n+1;
                if n == 1
                    lastj=j;
                end
                if n==40 && j>n+40 
                    lastj=j;
                    tst=timeev(j-80);  
                    starti=j-80;
                    1;
                    break
                end
                if n==40 && j<n+40 
                    lastj=j;
                    tst=timeev(j-39); 
                    starti=j-39;
                    2;
                    break  
                end
              
            end 
        end
       
        %(2)
        n=0;
        m=0;
        if tst==0 || starti >1500
            for j=1:length(irrev)
                if diff(j)>(2*irrstd)
                    n=n+1;
                    if n==30 && j>n+30
                        tst=timeev(j-60);
                        starti=j-60;
                        3;
                        break  
                    end
                    if n==30 && j<n+30
                        tst=timeev(j-29); 
                        starti=j-29;
                        4;
                        break
                    end
                end
            end
        end
        
        %(3)
        n=0;
        m=0;
        if tst==0 || starti >1500
            for j=1:length(irrev)
                if diff(j)>(2*irrstd)
                    n=n+1;
                    if n==20 && j>n+20
                        tst=timeev(j-40); 
                        starti=j-40;
                        5;
                        break
                    end
                    if n==20 && j<n+20
                        tst=timeev(j-19);  
                        starti=j-19;
                        6;
                        break
                    end
                end
            end
        end
        
        %(4)
        
        if tst==0 || starti >1500
            for j=1:length(irrev)
                if diff(j)>(2*irrstd)
                    n=n+1;
                    if n==10 && j>n+10
                        tst=timeev(j-20); 
                        starti=j-20;
                        7;
                        break
                    end
                    if n==10 && j<n+10
                        tst=timeev(j-9);  
                        starti=j-9;
                        8;
                        break
                    end
                end
            end
        end
        
    
    
        %at this point, j will be the index at which the criteria (whether
        %in 1, 2, or 3) was satisfied.  Carry this along: 
        startj = j;
        
        %(5)-(7) account for if the found start time was already counted as an 
        %event. If it has been, it will have been stored in "starttimes" below already.
        %Then start from startj (the point at which the criteria was met for
        %the first start time) and move forwards.  For (4) through (6), simply
        %repeat the steps above, but start after the end time of the last
        %flare (calculated below as endj).  (5) and (6) lower
        %"number-of-points" criteria in the same way as (2) and (3)
        %if starti > length(timeev) || endj > length(timeev)
        if starti > length(timeev) %testing this, not using endj from last flare
            'No classification of flare!';
            %print('line295')
        
        %if the time of start is actually below the first element in
        %timeev (the window of time being studied), the start time will lie
        %outside.  Skip these events for now, still plotting solar quiet and 
        %irradiance, but not start or end times - but will need to be fixed later.
        elseif tst<timeev(1) 
            
            'Start too early!';
            %print('line304')
            
        else
            tst=timeev(starti);
            %Now the fun bit.  Find irradiance peak for well-behaved data,
            %first determining the next time that a light curve returns below
            %the solar quiet for a certain number of consecutive points (50, in
            %this case)
    
            %N.B. this first piece of logic only works if the start time 
            %was found using the first condition! (and doesn't work great if 
            %the increase in light curve is so fast that you're quickly on the decay phase)
            m=0;
            for j=(starti+40):length(diff)
                if i == 2032 
                    if diff(j)<0
                        m=m+1;
                        if m==30
                            endj=j;
                            
                            break
                        end
                    end
                end                    
                if diff(j)<0
                    m=m+1;
                    if m==50
                        endj=j;
                        
                        break
                    end
                end
            end
            %testing this to account for no endj - add to allflares script?
            if endj>length(timeev)
                disp('no end found')
                endj=length(timeev)
            end

            if isnan(endj)
                disp('endj is NaN')
            else
            
            %find the tend, the end time for the event, by identifying the
            %position within timeev (the window of time for the event being
            %studied) at which endj lies
            tend=timeev(endj);
            
            %now the window of the flare is consolidated into "window,"
            %consisting of the irradiance values between startj and endj
            window=irrev(starti:endj);
            windowt=timeev(starti:endj);
            sqwindow = sqev(starti:endj);
            sq=sqev(1:starti);
                    
            %find the average irradiance values within a 15-point spread around
            %each point (greater than 7 positions from the start and less than
            %7 positions from the end, of course) - this will give a better
            %idea of the actual average at the peak
            avgmeans=NaN(length(window)-7,1);
            for k=8:length(window)-7
                avgmeans(k-7)=nanmean(window(k-7:k+7));
            end
            
            %find the maximum in this average of the window
            maxind1=find(avgmeans==max(avgmeans));
            
            %find where the maximum lies in the original window, calculated at
            %the beginning of studying the flare
            maxind=starti+maxind1;
            
            %find the time of the maximum
            maxt=timeev(maxind);
            end
        end
        
        %find if the start time of this event is before the end time of the
        %last or the flare is really after this calculated window
        if i>1 
            if tst < endtimes_corr(t-1,1) %|| max(irrev(endj:end))-max(window) > 2*irrstd
                
                %reset tst
                tst=0;
                
                %(5)
                for j=endj:length(irrev)
                    if diff(j)>(2*irrstd)
                        n=n+1;
                        if n==40 && j>n+40
                            tst=timeev(j-80);
                            starti=j-80;
                            break
                        end
                        if n==40 && j<n+40
                            tst=timeev(j-39); 
                            starti=j-39; 
                            break
                        end
                    else
                        n=0;
                    end
                end
                
                %(6)
                if tst==0
                    for j=endj:length(irrev)
                        if diff(j)>(2*irrstd)
                            n=n+1;
                            if n==30 && j>n+30
                                tst=timeev(j-60);
                                starti=j-70;
                                break
                            end
                            if n==30 && j<n+30
                                tst=timeev(j-29);
                                starti=j-29;
                                break
                            end
                        else
                            n=0;
                        end
                    end
                end
    
                %(7)
                if tst==0
                    for j=endj:length(irrev)
                        if diff(j)>(2*irrstd)
                            n=n+1;
                            if n==20 && j>n+20
                                tst=timeev(j-40);
                                starti=j-40;
                                break
                            end
                            if n==20 && j<n+20
                                tst=timeev(j-19); 
                                starti=j-19;
                                break
                            end
                        else
                            n=0;
                        end
                    end
                end
                
                %(8)
                if tst==0
                    for j=endj:length(irrev)
                        if diff(j)>(2*irrstd)
                            n=n+1;
                            if n==10 && j>n+10
                                tst=timeev(j-20);
                                starti=j-20; 
                                break
                            end
                            if n==10 && j<n+10
                                tst=timeev(j-9); 
                                starti=j-9;
                                break
                            end
                        else
                            n=0;
                        end
                    end
                end
                if tst==0
                    starti=100000;
                end
            end
        end
        
        %at this point, j will be the index at which the criteria (whether
        %in 1, 2, or 3) was satisfied.  Carry this along: 
        startj = j;
        
        %(5)-(7) account for if the found start time was already counted as an 
        %event. If it has been, it will have been stored in "starttimes" below already.
        %Then start from startj (the point at which the criteria was met for
        %the first start time) and move forwards.  For (4) through (6), simply
        %repeat the steps above, but start after the end time of the last
        %flare (calculated below as endj).  (5) and (6) lower
        %"number-of-points" criteria in the same way as (2) and (3)
    
        
        
        %if the time of start is actually below the first element in
        %timeev (the window of time being studied), the start time will lie
        %outside.  Skip these events for now, still plotting solar quiet and 
        %irradiance, but not start or end times - but will need to be fixed later.
        if starti > length(timeev) || endj > length(timeev) 
            'No classification of flare!';
       
        
        
        elseif tst<timeev(1) 
            
            'Start too early!';
          
            
            
        else
            tst=timeev(starti);
            %Now the fun bit.  Find irradiance peak for well-behaved data,
            %first determining the next time that a light curve returns below
            %the solar quiet for a certain number of consecutive points (50, in
            %this case)
    
            %N.B. this first piece of logic only works if the start time 
            %was found using the first condition! (and doesn't work great if 
            %the increase in light curve is so fast that you're quickly on the decay phase)
            m=0;
            for j=(starti+40):length(diff)

                if i == 2032
                    if diff(j)<0.5*irrstd
                    m=m+1;
                        if m==30
                            endj=j;
                            
                            break
                        end
                    end
                end      
                if diff(j)<0.5*irrstd
                    m=m+1;
                    if m==50
                        endj=j;
                        i;
                        
                        break
                    end
                end
            end
            if isnan(endj)
                
                disp('endj is NaN')
            else
            
            %find the tend, the end time for the event, by identifying the
            %position within timeev (the window of time for the event being
            %studied) at which endj lies
            tend=timeev(endj);
            
            %now the window of the flare is consolidated into "window,"
            %consisting of the irradiance values between startj and endj
            window=irrev(starti:endj);
            windowt=timeev(starti:endj);
            sqwindow = sqev(starti:endj);
            sq=sqev(1:starti);
                    
            %find the average irradiance values within a 15-point spread around
            %each point (greater than 7 positions from the start and less than
            %7 positions from the end, of course) - this will give a better
            %idea of the actual average at the peak
            avgmeans=NaN(length(window)-7,1);
            for k=8:length(window)-7
                avgmeans(k-7)=nanmean(window(k-7:k+7));
            end
            
            %find the maximum in this average of the window
            maxind1=find(avgmeans==max(avgmeans));
            
            %find where the maximum lies in the original window, calculated at
            %the beginning of studying the flare
            maxind=starti+maxind1;
            
            %find the time of the maximum
            maxt=timeev(maxind);
            if length(maxt)>1
                maxt=maxt(1);
            elseif length(maxt)<1
                maxt=(tend+tst)/2;
            end
            end
        end
        
    
        
        %You're done with finding critical times for that flare (well, assuming
        %it was nicely behaved enough for the code)! 
        
        %Store the times in arrays initialized above
        starttimes_corr(t,1)=tst;

        endtimes_corr(t,1)=tend;
        maxtimes_corr(t,1)=maxt;
        starttimes_corr(t,2)=i;
        endtimes_corr(t,2)=i;
        maxtimes_corr(t,2)=i;
        

    
    %%%%%END ORIGINAL MODEL-BUILDING%%%%%
    
        %calculate the impulsiveness value! Now the window is between the start
        %and end - repeat finding the avgmeans/average values with a 15 point
        %window - this can be changed to a smaller or larger window as you
        %wish.
        if endj < length(irrev) && starti < endj 
            window=irrev(starti:endj);
            windowt=timeev(starti:endj);
            sqwind=sqev(starti:endj);
            errwind=errev(starti:endj);
            sqpre=irrev(1:starti);
            sqerr = nanstd(sqpre);
        else

            continue
        end
        
        if length(windowt) < 12
            continue
        end

        windowthr=windowt(1):((windowt(2)-windowt(1))/8):windowt(length(windowt));
        
        avgmeans=NaN(length(window)-7,1);
        sqav = mean(sqev);
        for k=8:length(window)-7
            avgmeans(k-7)=nanmean(window(k-7:k+7));
        end
        models = NaN(length(windowt),1);
    
        %identify where the maximum values actually lie...
        maxI_304 = max(window);
        
        %...and find the index corresponding to this in the window...
        %(1) option 1 is to use the raw maximum of the dataset
        maxind1=find(window==maxI_304);
        maxind=starti+maxind1;
        %(2) option 2 is to use the average maximum over 15 point window
    %     maxind1=find(avgmeans==max(avgmeans));
    %     maxind=startj+8+maxind1;
        %...but also identify where that is in the original sqev, irrev, and
        %timeev arrays.
        
        sq_at_max = sqev(maxind);
        peakt304 = timeev(maxind);
        
        %now the irradiance at half height (relative to solar quiet) is just
        %the average of the solar quiet value at the max and the max irradiance
        %itself - this needs some work, consider other options.
        irr_hh = (sqav+maxI_304)/2;
        
        %if this half height actually has an irradiance value in it, find the
        %times corresponding to irradiance points which lie roughly around this
        %irradiance, and take the mean time in this window to find the low and
        %high values of time at half height
        if isempty(irr_hh)==0
            %base impulsiveness calculation only on the spline - see original
            %script for other options
            for k=2:length(window)
                if isnan(window(k))
                    window(k)=window(k-1);
                end
            end
            for k=length(window)-1:-1:1
                if isnan(window(k))
                    window(k)=window(k+1);
                end
            end
            spl = interp1(windowt,window,windowthr);
            splsq = interp1(windowt,sqwind,windowthr);
            splerr = interp1(windowt,errwind,windowthr);
            
            
            smspl = smooth(spl,5);
            vsmspl = smooth(spl,20);        
            
            %define error of each data point as the residual between a very
            %smoothed curve with the only slightly smoothed curve
            sqspl = smooth(splsq,5);
            errorbars = abs(vsmspl-spl');
            errm = nanmean(errorbars);
       	    for o=1:length(errorbars)
                if errorbars(o) == 0
                    errorbars(o) = errm;
                end
            end
            resid=smspl-sqspl;
            lenres=1:length(resid);
             %identify where the maximum values actually lie...
             
                maxI_304 = max(smspl);
                maxI_3042 = max(window);
    
                %...and find the index corresponding to this in the window...
                %use the raw maximum of the dataset - see other options in
                %original script
                maxind1=find(smspl==maxI_304);
                maxind2=find(window==maxI_3042);
                maxind=starti+maxind2;
    
            if maxind1 >1
               
                %...but also identify where that is in the original sqev, irrev, and
                %timeev arrays.
                
                sq_at_max = sqev(maxind);
                peakt304 = timeev(maxind);
    
                %now the irradiance at half height (relative to solar quiet) is just
                %the average of the solar quiet value at the max and the max irradiance
                %itself - this needs some work, consider other options.
                irr_hh = (sqav+maxI_304)/2;
    
                [~,ix] = min(abs(smspl(1:maxind1)-irr_hh));
                t_hh_low = windowthr(ix);
    
                [d,ix] = min(abs(smspl(maxind1:length(smspl))-irr_hh));
                t_hh_high = windowthr(maxind1+ix);
    
                if length(t_hh_low) > 1
                    t_hh_low = t_hh_low(1);
                end
    
                if length(t_hh_high) > 1
                    t_hh_high = t_hh_high(1);
                end        
                %distance between these two will be the FWHH
                t_hh = t_hh_high-t_hh_low;
    
                curly_I = maxI_304/t_hh;

                curly_I_2=(maxI_304-sq_at_max(1))/sq_at_max(1)/t_hh;
                
                curly_Is_best(t,1)=log(curly_I*1e6/24/3600);

                curly_Is_relative_best(t,1)=log(curly_I_2/24/60);

                curly_Is_best(t,2) = i;
                curly_Is_relative_best(t,2) = i;

                %redo modeling process with the high-res window to compare to
                %spline
    
    
                risehr_t = windowthr(1:maxind1)';
                decayhr_t = windowthr((maxind1+1):length(spl))';
                risehr = spl(1:maxind1)';
                decayhr = spl((maxind1+1):length(spl))';
                riseerr=splerr(1:maxind1)';
                decayerr=splerr((maxind1+1):length(splerr));
                risehr_mirror=NaN(2*length(risehr),1);
                risehr_mirrort=NaN(2*length(risehr_t),1);
                dt = risehr_t(2) - risehr_t(1);
    
            end
                
        end
        
        % clf
        % f2 = figure('visible','off');
        % plot(timeev,irrev)
        % hold on
        % xline(timeev(starti))
        % xline(timeev(endj))
        % saveas(gcf,'/Users/coletamburri/Desktop/bestflares/'+[string(i)+'.png'])

        event_curves(t,1:length(irrev))=irrev/(1e-3);
        event_times(t,1:length(timeev))=timeev;

         end

        
end


curly_Is_best(:,2) = sort(is);
curly_Is_relative_best(:,2) = sort(is);

    