function pathStr = projectPath(varargin)
%PROJECTPATH Build a path relative to the repository root.
pathStr = fullfile(projectRoot(), varargin{:});
end
