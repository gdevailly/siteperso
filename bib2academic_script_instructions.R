library(bib2academic) # remotes::install_github("petzi53/bib2academic", build_vignettes = TRUE)

# devailly.bib from Zotero (all document in collection "mes publications")
# select all -> right click -> export documents to bibtext
# TODO : rename Heat*seq token to "devailly_heatstarseq_2016"

bib2acad("devailly.bib", copybib = TRUE, abstract = TRUE)

# TODO : "If everything went smoothly make a copy of your static/files/citations/
# and content/publication folder.
# Then copy your .bib files into static/files/citations/
# and your .md files into content/publication."

# TODO : change selected to true in the .md of a few publications
# TODO : add pictures
