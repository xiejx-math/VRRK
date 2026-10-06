function files = download_suitesparse(dataDir,matrixNames)
%DOWNLOAD_SUITESPARSE Download requested SuiteSparse Matrix Market files.

if nargin<1 || isempty(dataDir)
    dataDir=fullfile('data','suitesparse');
end
if nargin<2 || isempty(matrixNames)
    matrixNames=["ash958","well1033","illc1033"];
end
if ~exist(dataDir,'dir')
    mkdir(dataDir);
end

catalogNames=["ash958","well1033","illc1033"];
catalogGroups=["HB","HB","HB"];
files=cell(numel(matrixNames),1);

for matrixIndex=1:numel(matrixNames)
    matrixName=string(matrixNames(matrixIndex));
    catalogIndex=find(catalogNames==matrixName,1);
    if isempty(catalogIndex)
        error('Unknown SuiteSparse matrix: %s.',matrixName);
    end

    matrixFolder=fullfile(dataDir,matrixName);
    files{matrixIndex}=fullfile(matrixFolder,matrixName+'.mtx');
    if exist(files{matrixIndex},'file')
        continue;
    end

    archive=fullfile(dataDir,matrixName+'.tar.gz');
    url=sprintf('https://sparse.tamu.edu/MM/%s/%s.tar.gz',...
        catalogGroups(catalogIndex),matrixName);
    fprintf('Downloading %s\n',url);
    websave(archive,url);
    untar(archive,dataDir);

    if ~exist(files{matrixIndex},'file')
        error('Downloaded archive does not contain %s.',files{matrixIndex});
    end
end
end

