function [settings] = readSettings(dataDir)
    setFile = fopen([dataDir 'settings.bin'],'r');
    
    settings{1} = char(fread(setFile,1,'char*1'));  %simtype
    settings{2} = fread(setFile,1,'int');               %ngp
    settings{3} = fread(setFile,1,'int');             %tsmax
    settings{4} = fread(setFile,1,'double');           %lref
    settings{5} = fread(setFile,1,'double');             %dt
    settings{6} = fread(setFile,1,'double');             %dx

    if settings{1}=='t'
        settings{7} = fread(setFile,1,'double');         %Re
        settings{8} = fread(setFile,1,'double');       %uref
    elseif settings{1}=='v'
        settings{7} = fread(setFile,1,'double');          %a
        settings{8} = fread(setFile,1,'double');     %lratio
        settings{9} = fread(setFile,1,'double');         %Re
        settings{10} = fread(setFile,1,'double');         %nu
    end

    fclose(setFile);
end
