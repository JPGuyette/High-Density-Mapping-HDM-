function hdm_processvid(dataset,ext,rotangle,crop)
    %% Process function inputs
        % dataset
        vidfiles=dir([dataset,'\*',ext]);
        src_dir=[dataset,'\'];
        if isempty(vidfiles)
            error('No videos found in the specified dataset')
        end
        vid=VideoReader([src_dir,vidfiles(1).name]);
        % Get video information
        numframes=vid.NumberOfFrames;
        framerate=vid.FrameRate;
        disp(['Frame Rate is:',num2str(framerate)]);
        disp(['Number of Frames is:',num2str(numframes)]);
        curim=read(vid,1);
        if size(curim,3)==3;
            curim=rgb2gray(curim);
        end
        % Image cropping
        if crop==1
            fig=figure;
            
            [curim,croprect1]=imcrop(curim);
            close(fig);
            croprect2=0;
        else
            croprect1=0;
            croprect2=0;
        end
        % Rotation angle    
        if isempty(rotangle)
            fig=figure;
            imshow(curim)
            h=imline;
            pos=wait(h);
            close(fig);
            dxy=pos(2,:)-pos(1,:);
            ang=atand(dxy(2)/dxy(1));
            if crop==1;
                curim=imrotate(curim,ang);
                fig=figure;
                [~,croprect2]=imcrop(curim);
                close(fig);
            end
        else
            ang=rotangle;
            if crop==1 && ang~=0;
                curim=imrotate(curim,ang);
                fig=figure;
                [~,croprect2]=imcrop(curim);
                close(fig);
            end
        end
        % Output folder
        out_dir=[src_dir,'images\'];
        mkdir(out_dir);
        
        % Write information to file
        save([out_dir,'imageinfo.mat'],'dataset','croprect1',...
            'croprect2','ang','numframes','framerate');
        
    %% Process and output images
        tic;
        starttime=now();
        bar=waitbar(0,'Processing video frames');
        parfor i = 1:numframes
            % Load image
            curim=read(vid,i);
            if size(curim,3)==3;
                curim=rgb2gray(curim);
            end
            % Crop image
            if crop==1
                curim=imcrop(curim,croprect1);
            end
            % Rotate image
            if ang~=0
                curim=imrotate(curim,ang);
            end
            % Crop image again
            if crop==1 && ang~=0
                curim=imcrop(curim,croprect2);
            end
            % Write image to file
            imwrite(curim,[out_dir,dataset,'_',numpad(i,4),'.tif'])
            t=now-starttime;
            x=i/numframes;
            timeremaining=datestr((t*(1-x))/x,'HH:MM:SS');
            waitbar(x,bar,['Approximate time remaining: ',timeremaining]);
            
        end
        close(bar);
        toc
end