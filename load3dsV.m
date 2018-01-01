function [X, Y, V, Z, data, header, par] = load3dsV(fn)

% --- load3dsV is an extend version of the original Nanonis routine to read
% the .3ds data files. it outputs X, Y meshgrids and V vector (this is the
% difference from the load3ds version) as well as data matrix which
% contains all recoreded channels for each X, Y and V point, header
% structure, and parameter matrix, as specified in the header. 
% load3ds  Nanonis 3ds file loader
%   [header, data, par] = load3ds(fn, pt_index) reads a Nanonis 
%   3ds file fn.
%   Without a second argument, only the header is returned.
%   With a second argument pt_index, the return value 'data'
%   contains the spectroscopy data set of point pt_index (zero
%   based). Each column in the data array corresponds to an
%   acquired channel.
%   The 'par' variable will contain the parameters of the data
%   set at the specified point.

data=''; header=''; par='';

if exist(fn, 'file')
    fid = fopen(fn, 'r', 'ieee-be');    % open with big-endian
else
    fprintf('File does not exist.\n');
    return;
end

% read header data
% The header consists of key-value pairs, separated by an equal sign,
% e.g. Grid dim="64 x 64". If the value contains spaces it is enclosed by
% double quotes (").
while 1
    s = strtrim(fgetl(fid));
    if strcmp(upper(s),':HEADER_END:')
        break
    end
    
    %s1 = strsplit(s,'=');  % not defined in Matlab
    s1 = strsplit_i(s,'=');

    s_key = strrep(lower(s1{1}), ' ', '_');
    s_val = strrep(s1{2}, '"', '');
    
    switch s_key
    
    % dimension:
    case 'grid_dim'
        s_vals = strsplit_i(s_val, 'x');
        header.grid_dim = [str2num(s_vals{1}), str2num(s_vals{2})];
        
    % grid settings
    case 'grid_settings'
        header.grid_settings = sscanf(s_val, '%f;%f;%f;%f;%f');
         
    % fixed parameters, experiment parameters, channels:
    case {'fixed_parameters', 'experiment_parameters', 'channels'}
        s_vals = strsplit_i(s_val, ';');
        header.(s_key) = s_vals;
        
    % number of parameters
    case '#_parameters_(4_byte)'
        header.num_parameters = str2num(s_val);
        
    % experiment size
    case 'experiment_size_(bytes)'
        header.experiment_size = str2num(s_val);

    % spectroscopy points
    case 'points'
        header.points = str2num(s_val);

    % delay before measuring
    case 'delay_before_measuring_(s)'
        header.delay_before_meas = str2num(s_val);
    
    % other parameters -> treat as strings
    otherwise
        s_key = regexprep(s_key, '[^a-z0-9_]', '_');
        header.(s_key) = s_val;
    end
end

% read the data 

    exp_size = header.experiment_size + header.num_parameters*4;
    %prepares data matrix data (X x Y x V x # of channels),
    %parameters matrix par (X x Y x # of parameters) and final X, Y (meshgrids) and V (vector)
    data=zeros(header.grid_dim(2),header.grid_dim(1),header.points,prod(size(header.channels)));
    par=zeros(header.grid_dim(2),header.grid_dim(1),header.num_parameters);
    Xmesh=zeros(header.grid_dim(2),header.grid_dim(1));
    Ymesh=zeros(header.grid_dim(2),header.grid_dim(1));
    sizeD=size(data);
    size_par=size(par);

    %get current position in file
    pos = ftell(fid);
          
    for j=0:header.grid_dim(2)-1
        for i=0:header.grid_dim(1)-1
            pt_index=j*header.grid_dim(1)+i;
            fseek(fid,pos,-1);
            fseek(fid, pt_index*exp_size, 0);

            temp1 = fread(fid, header.num_parameters, 'float');
            size_temp1=size(temp1);
            if size_par(3)==size_temp1(1)
                par(j+1,i+1,:)=temp1;
            end
            clear temp1
            
            temp2=fread(fid, [header.points prod(size(header.channels))], 'float');
            size_temp2=size(temp2);
            if sizeD(3)==size_temp2(1) &&  sizeD(4)==size_temp2(2)
                data(j+1,i+1,:,:) = temp2;
            end
            clear temp2             
        end
        
    end
    % Build X and Y vectors 
    X=linspace(0,header.grid_settings(3),header.grid_dim(1));
    Y=linspace(0,header.grid_settings(4),header.grid_dim(2));
    % Build the bias vector
    Vstart=par(1,1,1);
    Vend=par(1,1,2);
    % request user input on voltage divider ratio
    Vdivider = 1;  % input('Please specify voltage divider ratio \n');
    V=linspace(Vstart*Vdivider,Vend*Vdivider,header.points);
    
    % Build Z matrix
    Zind= find(ismember(header.experiment_parameters,'Z (m)')); %finds the channel of Z data
    Zind=Zind+length(header.fixed_parameters);  % corrects its index in par by adding the number of fixes parameters (see header structure)
    Z=par(:,:,Zind); 
    

fclose(fid);

end  % of function load3ds


function s = strsplit_i(str, del)
    s = {};
    while ~isempty(str),
        [t,str] = strtok(str, del);
        s{end+1} = t;
    end
end

 
