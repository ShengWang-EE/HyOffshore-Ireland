function root = projectRoot()
%PROJECTROOT Return the absolute path to the repository root.
root = fileparts(mfilename('fullpath'));
end
