

function mapout=jtfheatmap(mapname,size)
    if strcmp(mapname(1),'-')
        if strcmp(mapname(2),'-')
            mn=mapname(3:end);
            mapname=mapname(3:end);
        else
            mn=mapname(2:end);
        end
    else
        mn=mapname;
    end
    
    mapout=eval([mn,'(',num2str(size),');']);
    
    if strcmp(mapname(1),'-')
        mapout=flipud(mapout);
    end
end