function chkdir
err=0;
    if strcmp(pwd,'\\research.wpi.edu\gaudettelab\Projects')
       x=questdlg(['MatLab is currently working in the Gaudette Lab Projects folder.',...
           '  Typically analysis is performed in a subfolder of Projects.',...
           ' Would you like to continue?'],...
           'Screw up data management for the entire lab?',...
           'Continue ruining data management for the entire lab',...
           'Cancel and manually change the working directory',...
           'Cancel and manually change the working directory');
       if strcmp(x,'Continue ruining data management for the entire lab')
           x=questdlg(['Seriously?.',...
           '  You really want to F*** everything up for the whole lab.',...
           ' because your too !#%$@%#$^ lazy to change the directory?'],...
           'Screw up data management for the entire lab?',...
           'Yes, F*** it all',...
           'Fine, I''ll change directories',...
           'Fine, I''ll change directories');
            if strcmp(x,'Yes, F*** it all')
                compname = getenv('COMPUTERNAME');
                compname = lower(deblank(compname));
                % Get the Last name of the current user
                curuser=[getenv('userdomain'),'\',getenv('username')];
                % Determine target build ID
                if ~exist('buildin','var')
                    [buildin,buildyear]=GCE_curbuild;
                end
                buildid=[buildyear,'.',buildin]; 
               % Create session log file
                sessdate=datestr(now(),'yyyy_mm_dd_HH:MM:SS');
                %logs=dir([GCEloc,'~GCE\~log\',sessdate,'*']);
                fpath='f---info.txt';
                logfile=fopen(fpath,'at+');
                fprintf(logfile,'\n%s','############');
                fprintf(logfile,'\n%s',sessdate);
                fprintf(logfile,'\n%s',['Build: ',buildid]);
                fprintf(logfile,'\n%s',['User: ',curuser]);
                fprintf(logfile,'\n%s',['Workstation: ',compname]);
                fprintf(logfile,'\n%s',['Matlab Version: ',version]);
                fclose('all');
                x=0;
                bar=waitbar(0,'F***ing up the lab. Please wait...');
                while x<30
                    pause(1)
                    x=x+1;
                    waitbar(x/30,bar);
                end
            else
                err=1;
            end
       else
           err=1 ;
       end
       if err==1
          error('User decided not to F*** up the lab'); 
       end
    end
    clearvars compname curuser buildin buildyear sessdate x y
end