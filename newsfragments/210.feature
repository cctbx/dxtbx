FormatMultiImage: When constructing an imageset with the indices of some (not
all) single images in the container, we skip reading models for the images that
were not requested. In some cases this speeds up imageset construction by 8x.
