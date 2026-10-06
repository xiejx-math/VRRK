function A = mmread_simple(file)
%MMREAD_SIMPLE Read real coordinate Matrix Market matrices.
fid=fopen(file,'r');
assert(fid>=0,'Cannot open %s',file);
cleaner=onCleanup(@() fclose(fid));
header=lower(fgetl(fid));
tokens=strsplit(strtrim(header));
assert(numel(tokens)==5 && strcmp(tokens{1},'%%matrixmarket') && ...
    strcmp(tokens{2},'matrix') && strcmp(tokens{3},'coordinate') && ...
    ismember(tokens{4},{'real','integer','pattern'}) && ...
    ismember(tokens{5},{'general','symmetric'}),...
    'VRRK:MatrixMarketFormat','Only real/integer/pattern coordinate general or symmetric matrices are supported.');
symmetric=strcmp(tokens{5},'symmetric');
pattern=contains(header,'pattern');
line=fgetl(fid);
while startsWith(strtrim(line),'%')
    line=fgetl(fid);
end
d=sscanf(line,'%d %d %d');
m=d(1);
n=d(2);
nz=d(3);
if pattern
    C=textscan(fid,'%f %f',nz);
    v=ones(nz,1);
else
    C=textscan(fid,'%f %f %f',nz);
    v=C{3};
end
i=C{1};
j=C{2};
if symmetric
    off=i~=j;
    A=sparse([i;j(off)],[j;i(off)],[v;v(off)],m,n);
else
    A=sparse(i,j,v,m,n);
end
end
