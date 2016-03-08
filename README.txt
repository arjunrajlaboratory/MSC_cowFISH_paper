### README first

Welcome to the cowFISH repository!

FOLDER CONTENTS: 

all_csvs - metadata spreadsheet (contains all info for all extracted spot counts), all extracted spot counts from RNA FISH data (organized by experiment id, more info in metadata), all other data that can be contained in spreadsheet format (some cell diameter data, some manual tracking info for timelapse-FISH matchup)

data_import - contains scripts to import all of the data from all_csvs into R into large, metadata matched spreadsheets

extract_from_raw_data - matlab script that pulls spot counts from raw data and export to csvs (NOTE: no raw data contained in the repository, this is just to trace which spot counts came from which raw images)

FigureN - the final figure in PDF form, the illustrator file for each figure, and the scripts necessary to make the graphs or crop example images in each figure

Stats - supplemental statistics data generated in SAS

TO GENERATE GRAPHS FOR ALL FIGURES:

First, in R, set the working directory to the cowFISH_paper repository directory (using command setwd()).

Then, 

source(./master_figure_generator.R) 

This will generate graphs for each figure in a folder named the current date within the graphs subfolder of each figure folder.

TO GENERATE GRAPHS FOR AN INDIVIDUAL FIGURE:

First, in R, set the working directory to the cowFISH_paper repository directory (using command setwd()).

Then, 

source(./FigureN/FigureN.R) 

where 'N' is the figure number of interest. 

This will generate graphs for FigureN in a folder named the current date within the graphs subfolder of FigureN.

TO GENERATE CROPPED IMAGES FOR FIGURE 1:

Run the merge_and_crop_Fig1.m in MATLAB from within the ./Figure1/FISH_examples_images folder.

