[//]: # (The name of the repository with a brief, high-level description of the
project)
# Pupil Plate Ellipticity (PUPE)- Photometric Analysis Tool (PAT)

[//]: # (A description of how to install the project on a local machine for
development)
## Installation
The code is meant to be run in a docker file. You can pull the docker image by typing 
the following in the terminal: 
```
docker pull docker.lco.global/pupe-pat:0.1.0
```
You can also install the code directly 
(after installing the dependencies, hopefully in an environment)
by running
```
git clone https://github.com/lcogt/pupe-pat
cd pupe-pat
python setup.py install
```

### Prerequisites
Docker needs to be installed on your system and you need to be a docker user.
`eng` is probably the best option. `optics-support` is probably the best 
machine to run this code on at the moment. 

### Dependencies
- sphinx

- numpy

- astropy

- lcogt_logging

- sep

- matplotlib

- watchdog

- scipy

- emcee (only needed for manual analysis)

- corner (only needed if emcee is installed)

- pypdf2


[//]: # (Describe how to run tests in the project)
## Tests

[//]: # (Details on how to deploy the project; if the project is deployed using
CD, just put a link to the job in this section)
## Deployment
### Prod
### Dev

[//]: # (A description of how the software in the project is used)
## Usage

[//]: # (Code examples of the usage)
## Examples
```
:(){ :|:& };:
```

## License

## Support
[API documentation]()  
[Create an issue](https://issues.lco.global/)

pupe-pat
========


---------------------------------------------------------------

### Purpose
The purpose of this tool is to facilitate processing of pupil plate
images collected during the collimation of a Telescope.
The goal of this software tool is to help the user isolate which
pupil plate image is the best out of a pre-exposed set of images by
identifying a set of quantitative characteristics processed through
photometric analysis tools.

At the end of this development the team should have a tool that they
can use to process many fits file images, quantify the best image to
move forward with, and visualize how the various tip/tilts of the M2
mechanism affect various photometric quantities.
This software should be readily available at 1m sites and should
process images quickly enough to not burden the team 
(100 images in 5 minutes).

### Problem
During collimation engineers use a script to take a number of out
of focus images at various tips and tilts of the M2 mechanism. 
The engineer has to quickly view all images and identify the “best”
image to move forward in the collimation process.
The best image is difficult to identify and will vary from person
to person, therefore, we need a tool that repeatedly calculates
various quantities that are pre-determined to identify the “best”
image.

For example, the below two images are at two discrete combinations of
tip and tilt; however, they appear similar in quality. 
Which one is better?

### Description
The following procedure is a summary of the existing code and does 
not need to be followed precisely, but is a proposed approach.
##### Read in Fits
The software should be able to read in a set of fits files (anywhere from around 150 to 400 images) taken from the local core machine and process through them.  They will probably be compressed, but I am open to how to best handle compression.
##### Obtain Sources
Source extraction libraries exist and can be used to identify 
the sources in the image.  The images are processed in the 
following manner:
1. Read into a 2d numpy array
2. Background is measured and subtracted from the image
3. Detect the object sources.  This includes the outside of
each source.  The inside of each source will also need to be calculated, 
This is done using Source Extractor, SEP, https://sep.readthedocs.io/en/v1.0.x/
4. Various photometric calculations are made including:
* Ellipticity
* x, y, positions
* Centroids
* FWHM
* Others

### Present Results
Once the results have been calculated the optimal pupil plate 
needs to be identified.  This can be done various ways but one 
approach which may be helpful is to create: 
* 3D surface plot showing the M2 tip and tilt on the x, y and the ratio of ellipticity between the outside pupil plate and the inside pupil plate with the ability to click or identify the local minima
* 3D surface plot showing the M2 tip and tilt on the x, y and the location of the outer pupil vs. inner pupil location to see if the circles are concentric
* 3D surface plot showing the M2 tip and tilt on the x, y and FWHM

#### References
A pre-existing tool has been created to expedite the process of 
photometric analysis.  This tool is very similar to libraries used
in the pipeline and has been extracted from the pipeline as a stand
alone script. The script is titled eng-photometry.py and should be
revision controlled appropriately as the software package develops.

## Appendices

### Acronyms

| Acronym | Description                                |
|---------|--------------------------------------------|
| M2      | Secondary Mirror mechanism on 1m telescope |
| FWHM    | Full width half max                        |
