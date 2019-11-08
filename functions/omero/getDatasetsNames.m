function datasetsNames = getDatasetsNames(datasetsLinks)
	datasetsNames = {};
	for j=0:datasetsLinks.size()-1,
        dataLink = datasetsLinks.get(j);
        dataset = dataLink.getChild();
        datasetName = strtrim(char(dataset.getName().getValue()));
        if (isempty(datasetName))
            datasetName = 'no name';
        end
        datasetsNames{end+1} = datasetName; 
    end
end
