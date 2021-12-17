# EDIBLES_Spectrum -- Haoyu Version
This is an integrated tool box for the EDIBLES survey. This project takes place at European Southern Observatory (ESO) in Chile and uses the 8-m class Very Large Telescope (VLT) that are the world's leading optical telescope. The EDIBLES survey acquired high-resolution spectra of 123 stars, and is aimed to solve the long-standing mystery in interstellar medium known as the Diffuse Interstellar Bands (DIBs).

This toolbox is to select, process, presents, and analysis the spectral data of the EDIBLES project. In this repo is the branch maintained by Haoyu Fan who worked with the EDIBLES project as a postdoc between 2020 and 2021. While many of the functions are also available in the main branch, the edition is more about uniforming the inputs and outputs of different functions so they can talk to each other, and thus "integrated". Haoyu also added some new function to allow point-accuracy and interactive processing of the data, so all the tools can be used in console or in script. 

Some of the major functions include:

- Select fits file based on target, wavelength, observation time etc.
- Read in single or multiple fits file and store the data cubic.
- Shift, align, trim the data cube. Convert the units of wave_grid and flux.
- Filter and mask bad data points.
- Normalize the with multiple methods and models.
- Build customized data of the absorption lines and fit to the data.
- Presents the spectral data with great flexibility, e.g. with fitted model, add highlighted species, user defined x/y labels.
- De-noise and output the spectral data as array (for further analysis).
- Advanced models for specific projects.