function imagesNames = getImagesNames(imagesLinks)
	imagesNames = {};
    for j = 0:imagesLinks.size()-1,
        imagesNames{end+1} = strtrim(char(imagesLinks.get(j).getName().getValue()));
    end
end
    