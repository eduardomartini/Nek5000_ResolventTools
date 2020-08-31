params = dlmread("params.input");
tol    = params(6);

files = strsplit(fileread(fileList));
files = files(1:end-1);

timeNorm=load(normFile);
timeNorm = timeNorm( timeNorm(:,3)==1,:);  %keep only timesteps writen to disk
p = find(timeNorm(:,2)>tol*max(timeNorm(:,2)));



if isempty(p)
      exit  
else
      p = p([1,end]);
end

disp(['Writting files to rename to "filesToRename.txt"']);
fid = fopen('filesToRename.txt','w');
for ifile = [1:p(1)-1,(p(2)+1):size(timeNorm(:,2))]
    disp(['    ' files{ifile}]);
    fprintf(fid, '%s\n', files{ifile});
end
fclose(fid)


