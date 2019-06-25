function [t,positions] = read_molPos3(filename,nframes)

positions.Ga = cell(nframes,1);
positions.Gi = cell(nframes,1);
positions.Ra = cell(nframes,1);
positions.Ri = cell(nframes,1);
t=nan(nframes,1);
fid=fopen(filename);

frameid = 1;
while ~feof(fid)
    % each line contains:
    % t x1 y1 x2 y2 .. xn yn
    % for all n species of the particular time.
    % all 4 species are listed as Ga Gi Ra Ri
       
    currline = fgetl(fid);
    [currt,x,y,z]=entry_to_xyz(currline);
    if ~isnan(currt) && isnan(t(frameid))
        t(frameid)=currt;
    end
    positions.Ga{frameid} = [x,y,z];
    
    currline = fgetl(fid);
    [currt,x,y,z]=entry_to_xyz(currline);
    if ~isnan(currt) && isnan(t(frameid))
        t(frameid)=currt;
    end
    positions.Gi{frameid} = [x,y,z];
        
    currline = fgetl(fid);
    [currt,x,y,z]=entry_to_xyz(currline);    
    if ~isnan(currt) && isnan(t(frameid))
        t(frameid)=currt;
    end
    positions.Ra{frameid} = [x,y,z];
    
    currline = fgetl(fid);
    [currt,x,y,z]=entry_to_xyz(currline);
    if ~isnan(currt) && isnan(t(frameid))
        t(frameid)=currt;
    end
    positions.Ri{frameid} = [x,y,z];
    frameid=frameid+1;
end
fclose(fid);

end

function [t,x,y,z]=entry_to_xyz(line)
    coords=sscanf(line,'%g'); 
    if numel(coords)>1
        t = coords(1);
        x = coords(2:3:end);% skip the time value, which is the first entry
        y = coords(3:3:end);
        z = coords(4:3:end);
    else
        t = nan; x=[]; y=[]; z=[];
    end
end