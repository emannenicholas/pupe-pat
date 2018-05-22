[//]: # (The name of the repository with a brief, high-level description of the
project)
# Pupil Plate Ellipticity (PUPE)- Photometric Analysis Tool (PAT)

[//]: # (A description of how to install the project on a local machine for
development)

## Purpose
The purpose of this tool is to facilitate processing pupil plate
images collected during the collimation of a telescope. A "pupil plate"
is a defocused star that is defocused so much that it appears as a donut.
The shape of that donut can be used to characterize (and hopefully eliminate)
optical distortion.


The goal of this software tool is to help the user isolate which
pupil plate image is the best out of a pre-exposed set of images by
identifying the best position for the tip and tilt of the M2 mechanism.

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
Currently the tests are run manually. Our test data set was taken on
kb05 at bpl on 20180409. We use the whole raw directory of images.
To reduce the test data run the following on the `optics-support` machine
as the `eng` user:
```
docker run --rm -it -v /home/eng/pupe-pat-test-data:/input -v /home/eng/pupe-pat-tests-{0.1.0}:/output docker.lco.global/pupe-pat:{0.1.0} run_pupepat --output-dir /output --input-dir /input
```
Replace `{0.1.0}` with the current version.
The tests pass if all of the blue circles should align with the donuts
in the PDF files produced in the output directory.

[//]: # (Details on how to deploy the project; if the project is deployed using
CD, just put a link to the job in this section)
## Deployment
Currently new versions are built manually. You can run
```
docker build -t docker.lco.global/pupe-pat:{0.0.1} .
docker push docker.lco.global/pupe-pat:{0.0.1}
```
replacing `{0.0.1}` with the current version.

[//]: # (A description of how the software in the project is used)
## Usage
During collimation, engineers use a script to take a number of out
of focus images at various tips and tilts of the M2 mechanism. 
Previously, the best image was identified by eye.
This code is meant to quantify this project.


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

We currently only support the LCO 1-m class of telescopes.
[//]: # (Code examples of the usage)
## Examples
```
:(){ :|:& };:
```

### Acronyms

| Acronym | Description                                |
|---------|--------------------------------------------|
| M2      | Secondary Mirror mechanism on 1m telescope |
| FWHM    | Full width half max                        |

## License
This code is licensed under GPL v3.0. See License.md for more information.
## Support
 [Create an issue](https://issues.lco.global/)

---------------------------------------------------------------
