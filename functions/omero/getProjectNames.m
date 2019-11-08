function projectNames = getProjectNames(projects)
	projectNames = {};
	for i=0:projects.size()-1,
	    projectName = strtrim(char(projects.get(i).getName().getValue()));
        if (isempty(projectName))
            projectName = 'no name';
        end
        projectNames{end+1} = projectName;
    end
end
