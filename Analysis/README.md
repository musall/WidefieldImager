# Widefield analysis
Clone or download the repository to your analysis PC and add it to your Matlab path. You can analyze your own imaging data (ideally captured with [WidefieldImager](https://github.com/musall/WidefieldImager) or download a demo recording from [here](https://drive.google.com/drive/folders/1OUmZP7cyI6O6wAWcOY895BurBda8JANN).

There are two basic tutorials for basic analysis of imaging data (Tutorial_basicAnalysis.m) and dimensionality reduction (Tutorial_dimReduction.m). In either tutorial, first navigate to the folder that contains your imaging data. For demo recordings, first extract the zip file and then assign the path to the data as the variable 'dataPath'.

```dataPath = path to your data;```

The first tutorial 'Tutorial_basicAnalysis' will work on any of the example datasets but some of the basic variables (such as the name of the imaging files) need to be adjusted for different datasets. By default, the settings match the dataset in the preprocessed folder (preproc_tactile_hindpawMap). This dataset is the smallest and quickest to use to get an overview on how to analyze imaging data. When using any of the datasets in the 'raw-datasets' folder, make sure to adjust the name of the imaging files to 'Frames_2_640_540_uint16' and set opts.preProc to false. To clearly see visual responses in phasemaping datasets, you should also increase the duration of opts.postStim to 5 seconds or longer.

The second tutorial 'Tutorial_dimReduction' will work for any of examples in the 'raw-datasets' folder.

Lastly, computePhaseMapsRaw can be used to create visual phase maps for the localization of visual areas in cortex. Use either example dataset 2 or 3 from the 'raw-datasets' folder to test this function.

# Tutorial_basicAnalysis
This is a demo script for some basic analysis that can be done with widefield imaging data. It shows how to load imaging data and collect it in a larger data array for subsequent analysis. When using 'preproc_tactile_hindpawMap', imaging data is already properly aligned to a stimulus onset and contains single-channel imaging data. When using other example datasets that contain imaging data with blue and violet illumination, the data is split into blue and violet imaging stacks (in this case opts.preProc should be set to false). If opts.hemoCorrect is true, blue and violet data can also be used for hemo-dynamic correction of each pixel in the imaging data.

After loading the imaging data, it is spatially downsampled to reduce the data size and reduce the amount of single-pixel noise.

All example data sets contain a stimulus trigger and the script shows how to produce a trial-averaged stimulus response map and an averaged flourence trace from an area of interest.
Lastly, the script loads the trial-averaged data into a browsing GUI called 'compareMovie'. Here you can explore the imaging dataset by using the bottom slider to move through imaging frames. The GUI also contains several other useful options such as adjust the colormap and color range, rotating the image and isolating a flourescence trace (by clicking get trace and then clicking the part of the image to compute the trace in). It can also be used to save the imaging data as a movie file in .avi format or isolate different areas through tresholding.


![picture](https://github.com/musall/WidefieldImager/blob/master/images/compareMovie_example.png)

# Tutorial_dimReduction
This is a demo script to demonstrate a pre-processing pipeling for dimensionality reduction and subsequent hemodynamic reduction of widefield imaging data. 
Type 

```edit Tutorial_dimReduction```

to open the script. All variables and analysis steps are explained the file and should be largely self-explanatory. Please write an issue report if there are any problems and I'll try to assist.

The pipeline has multiple important steps: blockSVD separates the imaging frames into smaller blocks and performs linear dimensionality using randomized SVD. This steps returns block-wise data 'bV' and 'bU'.
A second SVD is then used to isolate common temporal dimensions across all blocks, representing the main temporal components for the whole data set. These components 'nV' are then used to compute the corresponding whole-frame spatial components 'U'.

Lastly, 'SvdHemoCorrect' performs hemodynamic correction on the low-dimensional data by regressing out fluoresence with violet illumination from blue illumination frames.

The resulting datafile 'Vc' represents the corrected imaging data. To restore individual frames simply transpose U and Vc (second dimension of Vc are frames). 
For example

```rawData = svdFrameReconstruct(U,Vc(:,1:10));```

will restore the first 10 frames in the current dataset.

The 'index frameCnt' contains the number of frames per trial and 'stimTime' indicates at which trial a stimulus was presented'. Using these variables you can reconstruct imaging data from different trials or responses to stimulus events of interest.

# computePhaseMapsRaw
This function generates visual phase maps from cortical imaging data. Download raw example dataset 2 or 3 or use your own cortical imaging data. A detailed description on visual phase maps and how to create the visual stimulation for visual mapping can be found in this [protocol](https://www.nature.com/articles/nprot.2016.158) by Juavinett et al.

To use the function, simply provide the path to the dataset ('dataPath') and the name of the animal (Plex66 for dataset2 or Fez73 for dataset3).

![picture](https://github.com/musall/WidefieldImager/blob/master/images/phaseMap_example.png)
