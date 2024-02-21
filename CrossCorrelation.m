% Mohit Bokariya
% Email : mohit.seismology@gmail.com
%After doing normalization and spectral whitening using norm_white1.sh 
%Cross Correlation and Stacking between two station of SAC data with multiple stations
%IN this code we will compute stack(GREEN Function) of all possible station pair
clear all; 
%ALL Station data folder List
stnlist=["H0010","H0100","H0200","H0290","H0400","H0500","H0600","H0700","H0800","H1110","H1170","H1240","H1330","H1405","H1470","H1540","H1580","H1630"];

%LOOP for all possible station pair

for s1=1:length(stnlist)-1;
     s2=stnlist(s1);
     for s3=s1+1:length(stnlist);
        s4=stnlist(s3);

        clearvars -except s1 s2 s3 s4 stnlist
        
        %Path for First station DATA file
        Path1 = '/media/iiser/Shaktimaan/Kathmandu/XF';
        cd(Path1);
        cd(s2);
        
        %Reading All SAC files and storing in DATA1 from Path1/s2
        DATA1 = dir('*.SAC');
        len1=length(DATA1);     %number of First station SAC files in given folder
        
        
        for i=1:len1
            fileNames1 = (DATA1(i).name);    %storing all file name one by one 
            fileNames1(end-11:end);             %storing Date part from file name
            f1=convertCharsToStrings(fileNames1(end-11:end));   %converting char to string
            date1(i)=append(f1);  %Storing All Dates of DATA1 first station files 
        
        end
        
        
        %Path for second station DATA file
        Path2 = '/media/iiser/Shaktimaan/Kathmandu/XF';
        cd(Path2)
        cd(s4)
        %Reading All second station SAC files and storing in DATA2 from Path2/s2
        DATA2 = dir('*.SAC');
        len2=length(DATA2);     %number of Second station SAC files in given folder
        for j=1:len2
            fileNames2 = (DATA2(j).name);    %storing all file name one by one 
            fileNames2(end-11:end);             %storing Date part from file name
            f2=convertCharsToStrings(fileNames2(end-11:end));       %converting char to string
            date2(j)=append(f2);   %Storing All Dates of DATA2 second station files
        
        end
        
        %Coommon Dates of Two Station Data sets
        CommonDates=intersect(date1,date2);
        len3=length(CommonDates);
        
        %IF there are no matching date between two stations then those station will be skip
        [r c]=size(CommonDates);
            if c==0
           
                continue
            end
            
        
        %Computing Common Field numbers of first station
        CommonFieldnum1=[];
        for i=1:len1
            fileNames1 = (DATA1(i).name);
            fileNames1(end-11:end);
            f1=convertCharsToStrings(fileNames1(end-11:end));
           
             for j=1:len3
        
                if f1==CommonDates(j)
                   
                  CommonFieldnum1=[CommonFieldnum1 i];
                end
             end
            
        end
        
        
        %Uncommon Field numbers which are not in Second station Data
        UncommonFieldnum1=setxor(CommonFieldnum1,1:len1);
        DATA1(UncommonFieldnum1)=[];
        
        %Computing Common Field numbers of Second station
        CommonFieldnum2=[];
        for i=1:len2
            fileNames2 = (DATA2(i).name);
            fileNames2(end-11:end);
            f2=convertCharsToStrings(fileNames2(end-11:end));
           
             for j=1:len3
        
                if f2==CommonDates(j)
                   
                  CommonFieldnum2=[CommonFieldnum2 i];
                end
             end
            
        end
        
        %%Uncommon Field numbers which are not in first station Data
        UncommonFieldnum2=setxor(CommonFieldnum2,1:len2);
        DATA2(UncommonFieldnum2)=[];
        
        
        %Computing Cross Correlation--
        Cross_Correlation=[];
        
        for i=1:len3
        
            cd(Path1)
            cd(s2)
            counts1=readsac(DATA1(i).name);
            
            cd(Path2)
            cd(s4)
            counts2=readsac(DATA2(i).name);
            
            [Cross_Correlation(i,:),lag]=(xcorr(counts1.DATA1,counts2.DATA1,300));
        end
        
        %Stacking Of Both station's Cross Correlations
        stack=0;
        for i=1:len3
            stack=stack+(Cross_Correlation(i,:)/len3);
        end

        %average of positive and negative side of all cross-correlations
        len4=length(stack);
        ln=(len4-1)/2;
        STACK=[stack(ln+1) flip(stack(1:ln))+stack(ln+2:len4)];
        STACK=STACK/2;  

        %Applying Bandpass on stack data
        
        bp=bandpass(stack,[0.0125 0.1],1);  %0.02 0.1 passband frequency range of the filter in hertz and 1 is frequency sample
        BP=bandpass(STACK,[0.0125 0.1],1);
        final_stack=bp/max(bp);   %normalization
        Final_stack=BP/max(BP);                  
        plot(lag,final_stack)

        %xlim([-1000 1000])
        
        
        %Giving name of final stack Result
        F1=fileNames1(1:end-12);  %first station- network.station.cmpound
        F2=fileNames2(1:end-13);   %second station- network.station.cmpound
        %Combining both station name
        filename=append(F1,F2,'.SAC');  % +ve and -ve side stack with .SAC extension
        Filename=append(F1,F2,'.SAC');    % average stack of +ve and -ve with .sac extension 
        
        stn1=F1(1:end-4);           %First network.station name
        stn2=F2(1:end-4);           %second network.station name
        stn=append(stn1,stn2);
        STN=append(stn(4:7),stn(11:13));
        

        cmp=F2(end-2:end);
        
        %Computin Azimuth angel and Distance between stations
        
        lat1=counts1.STLA;  %Latitude of first station
        lat2=counts2.STLA;   %Latitude of second station
        
        long1=counts1.STLO;   %Longitude of first station
        long2=counts2.STLO;    %Longitude of second station
        
        Azi=azimuth(lat2,long2,lat1,long1,'degree');        %Azimuth angle between two station in degree
        BackAzi=azimuth(lat1,long1,lat2,long2,'degree');     %Back Azimuth angle between two station in degree
        Dist=distance(lat1,long1,lat2,long2);    %distance between two station in degree
        Dist=deg2km(Dist);                        %distance between two station in KM
        
        
        %Path for storing Stack Data
        cd '/media/iiser/Shaktimaan/Kathmandu/Stack';
        
        
        %Saving Stack Data with headers
        %+ve and -ve stack data

        
        %Creating folder for normal stack with +ve and -ve side data
        if not(isfolder('Positive_Stack'))
            mkdir('Positive_Stack')
        end
        
        %Creating folder for average of both -ve and +ve side of cross-correlation
        if not(isfolder('Average_Stack'))
            mkdir('Average_Stack')
        end
        
        cd '/media/iiser/Shaktimaan/Kathmandu/Stack/Positive_Stack'
        mksac(filename,final_stack,lag(1),'EVLA',counts1.STLA,'EVLO',counts1.STLO,'STLA',counts2.STLA,'STLO',counts2.STLO,'DIST',Dist,'AZ',Azi,'BAZ',BackAzi,'KSTNM',STN,'KCMPNM',cmp)
        
        
        cd '/media/iiser/Shaktimaan/Kathmandu/Stack/Average_Stack'
        % average stack data
        mksac(Filename,Final_stack,lag(1),'EVLA',counts1.STLA,'EVLO',counts1.STLO,'STLA',counts2.STLA,'STLO',counts2.STLO,'DIST',Dist,'AZ',Azi,'BAZ',BackAzi,'KSTNM',STN,'KCMPNM',cmp)
        %disp(Filename)
        
    end
end



cd '/home/iiser/Batman/Codes'

