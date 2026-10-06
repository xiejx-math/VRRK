function folder = new_result_folder(parent,caseLabel)
%NEW_RESULT_FOLDER Preserve all prior runs, including repeated configurations.
if ~exist(parent,'dir')
    mkdir(parent);
end
[~,uniquePart]=fileparts(tempname(parent));
folder=fullfile(parent,char(string(caseLabel)+"_"+string(uniquePart)));
[ok,message]=mkdir(folder);
assert(ok,'VRRK:OutputDirectory','Cannot create output folder: %s',message);
end
