function archive = package_release()
%PACKAGE_RELEASE Build a source-only ZIP without local results or datasets.
root=fileparts(mfilename('fullpath'));
archive=fullfile(fileparts(root),'VRRK-codes-GitHub-release.zip');
% Explicit allowlist prevents results, downloaded data, and user paths leaking.
entries={'README.md','CHANGELOG.md','LICENSE','CITATION.cff',...
    '.gitignore','.gitattributes','data/README.md'};
patterns={'*.m','src/*.m','tests/*.m','tests/*.mtx'};
for patternIndex=1:numel(patterns)
    files=dir(fullfile(root,patterns{patternIndex}));
    for fileIndex=1:numel(files)
        absolute=fullfile(files(fileIndex).folder,files(fileIndex).name);
        entries{end+1}=absolute(numel(root)+2:end); %#ok<AGROW>
    end
end
zip(archive,entries,root);
fprintf('Source-only release: %s\n',archive);
end
