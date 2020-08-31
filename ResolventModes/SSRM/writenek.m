function status = writenek(fname,data,lr1,elmap,time,istep,fields,emode,wdsz,etag)
%
% This function writes binary data in the nek5000 new file format
%
%   status = writenek(fname,data,lr1,elmap,time,istep,fields,emode,wdsz,etag)
%
%   INPUT
%   - fname:  name of the file
%   - data:   nek5000 data ordered as (iel,inode,[x|y|(z)|u|v|(w)|p|T|s_i])
%   - lr1:    element-size vector (lx1,ly1,lz1)
%   - elmap:  reading/writing map of the elements in the file
%   - time:   simulation time
%   - istep:  simulation step
%   - fields: fields saved in the file
%   - emode:  endian mode 'le' = little-endian, 'be' = big-endian
%   - wdsz:   single (4) or double (8) precision
%   - etag:   tag for endian indentification
%
%   OUTPUT
%   - status: status (< 0 something went wrong)
%
%
% Last edit: 20151028 Nicolo Fabbiane (nicolo@mech.kth.se)
%

%--------------------------------------------------------------------------
% INITIALIZATION AND DATA CONSISTENCY
%--------------------------------------------------------------------------
%
% number of element, points-per-element, fields
[nel,npel,~] = size(data);
%
% check element dimension
if (npel ~= prod(lr1)), disp('npel ~= prod(lr1)'), status = -3; return, end
%
% check number of elements (TODO: multiple files not supported)
if (nel ~= size(elmap,2)), disp('nel ~= size(elmap,2)'), status = -4; return, end
%
% compute number of active dimensions
ndim = 2 + (lr1(3)>1);
%
% getfields [XUPT]
var=zeros(1,5);
if sum(fields == 'X') > 0
  var(1) = ndim;
end
if sum(fields == 'U') > 0
  var(2) = ndim;
end
if sum(fields == 'P') > 0
  var(3) = 1;
end
if sum(fields == 'T') > 0
  var(4) = 1;
end
if sum(fields == 'S') > 0
  var(5) = 0; % TODO: scalars not implemented
end
nfields = sum(var);
%
% check fields
if (size(data,3) ~= nfields), disp('size(data,3) ~= nfields'), status = -5; return, end
%
% set file id  
fid = 0; % TODO: multiple files not supported
%
% set tot number of files
nf = 1; % TODO: multiple files not supported

%--------------------------------------------------------------------------
%  OPEN FILE
%--------------------------------------------------------------------------
[outfile,message] = fopen(fname,'w+',['ieee-' emode]);
if outfile == -1, disp(message), status = -1; return, end

%--------------------------------------------------------------------------
% WRITE HEADER
%--------------------------------------------------------------------------
%
% write header
header = sprintf('#std %1i %2i %2i %2i %10i %10i %20.13E %9i %6i %6i %s\n',...
                 wdsz,lr1(1),lr1(2),lr1(3),nel,nel,time,istep,fid,nf,fields);
header(end+1:132) = ' ';
fwrite(outfile,header,'char');
%
% write endian tag
fwrite(outfile,etag,'float32');
%
% write element map
fwrite(outfile,elmap,'int32');

%--------------------------------------------------------------------------
% WRITE DATA
%--------------------------------------------------------------------------
%
% word size (double/single precision)
if (wdsz == 4)
    realtype = 'float32';
elseif (wdsz == 8)
    realtype = 'float64';
else
    fprintf('ERROR: could not interpret real type (wdsz = %i)',wdsz);
    status = -2; return
end
%
% write data
for ivar = 1:length(var)
    idim0 = sum(var(1:ivar-1));
    for iel = elmap
        for idim = (1:var(ivar))+idim0
            fwrite(outfile,data(iel,:,idim),realtype);
        end
    end
end

%--------------------------------------------------------------------------
% CLOSE FILE
%--------------------------------------------------------------------------
status = fclose(outfile);