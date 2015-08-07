function [myquery] = myinsitu2fn(filename, to_emesh_inp)
% [myquery] = myinsitu2fn(filename, to_emesh_inp)
% Converts filenames between insitu MySQL database and emesh
% to_emesh_inp can be omitted
% If to_emesh_inp is given and not equal 0, filename will be converted to
% enmesh compatible filename
% Otherwise filename will be converted from emesh to MySQL query compatible
% format
%

    if verLessThan('matlab', '7.2')
        error('This funciton uses ''regexpi'', which wasn''t implemented until Matlab 7.2');
    end

    to_emesh = 0;
    
    if exist('to_emesh_inp')
        if ~isempty(to_emesh_inp) || to_emesh_inp ~= 0,
           to_emesh = 1;
        end
    end
    
    if to_emesh == 0,
        % Convert filename to insitu MySQL compatible query 
        if strncmp(filename, 'insitu',6),
            fnum = regexpi(filename, '\d+', 'match');
            flen = length(fnum{1}) - 3;
            fdir = fnum{1}(1:flen);
            fdig = fnum{1}(flen+1:end);
            if str2num(fdig) == 0,
                fdir = str2num(fdir) - 1;
                fdir = int2str(fdir);
            end
        else
            fdir = '0';
        end
        
        fullname = strcat('img_dir_',fdir,'/',filename);
        myquery = fullname;
    else
        % Convert filename to emesh compatible
        fn_idx = regexpi(filename, '\/');
        myquery = filename(fn_idx+1:end);
    end

    
end
