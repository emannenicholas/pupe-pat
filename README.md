[//]: # (The name of the repository with a brief, high-level description of the
project)
# Pupil Plate Ellipticity (PUPE)- Photometric Analysis Tool (PAT)

[//]: # (A description of how to install the project on a local machine for
development)

## Purpose
The purpose of this tool is to facilitate processing pupil plate
images collected during the collimation of a telescope. A "pupil plate"
is a defocused star that is defocused so much that it appears as a donut.
The shape of that donut can be used to characterize
(and hopefully be able to eliminate) optical distortion.

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
docker pull docker.lco.global/pupe-pat:{0.1.0}
docker run --rm -it -v /home/eng/pupe-pat-test-data:/input -v /home/eng/pupe-pat-tests-{0.1.0}:/output docker.lco.global/pupe-pat:{0.1.0} run_pupepat --output-dir /output --input-dir /input
```
Replace `{0.1.0}` with the current version.
The tests pass if all of the blue circles should align with the donuts
in the PDF files produced in the output directory.

Eventually this procedure will be replaced with automatic testing in Jenkins.

[//]: # (Details on how to deploy the project; if the project is deployed using
CD, just put a link to the job in this section)
## Deployment
Currently new versions are built manually. You can run
```
docker build -t docker.lco.global/pupe-pat:{0.0.1} .
docker push docker.lco.global/pupe-pat:{0.0.1}
```
replacing `{0.0.1}` with the current version. Eventually we will replace
this procedure with continuous deployment in Jenkins.

[//]: # (A description of how the software in the project is used)
## Usage
During collimation, engineers use a script to take a number of out
of focus images at various tips and tilts of the M2 mechanism. 
Previously, the best image was previously identified by eye.
This code is meant to quantify this process. We currently only support
the LCO 1-m class of telescopes due to scaling the initial guess of the
donut size, but this should be straightforward to extend.

Images to be analyzed should be taken significantly defocused.
The stars should look like a donut with a hole in center. Currently,
we require that the inner hole have a radius of at least 5 pixels. 
The outer edge of the donuts should have a radius of at least 50 pixels
at the moment. This corresponds roughly to a defocus of +-12 mm for the 
Sinistro cameras. There are plans underway to shrink this requirement to a 
25 pixel radius which would allow defocus parameters down to +-6 mm.
This is currently a requirement in pixels rather than on defocus, so for the
sbigs in 1x1 binning we can currently get down to +-7.5 mm. 

The pupil plate images need to be taken on a bright, isolated star. A cut of 
V~4 seems to work well. This produces images that are high signal-to-noise
(half the saturation value) in about 1 second exposures. Fainter stars can be
used if they are isolated (the donuts do not overlap) and the exposure times
are long enough to get ~20,000 counts in the images.

The code automatically detects sources using Source Extractor in Python,
SEP, https://sep.readthedocs.io/en/v1.0.x. The code then makes a cutout
around the donut and fits a model using astropy.modeling that has one 
value in the "hole" of the  donut, one value for the value in the 
pupil plate and a background level. The centers and radii of the circles 
that form the boundaries of the pupil plate are also allowed to vary. 
There is code in this repo to fit ellipses instead of circles, and to fit
brightness gradient across the donut, but they are currently disabled.

The results are saved in the output directory specified by the user at runtime.
Each star in each image that is fit as a pupil plate produces a pdf file
with the filename of the image followed by a source id. At the end of the
analysis the pdfs are merged by defocus position into
pupe-pat-{focdmd}-donut.pdf (default) where {focdmd} is replaced 
with the requested focus of the images. The code also produces a "quiver"
plot (replacing donut with quiver in the name above). This shows the magnitude
and direction of the offset of the centers of the inner and outer edges of
the pupil plate. The "best" M2 tip and tilt are marked with a red X.

The code can be run in two modes (see below): real-time and offline modes.
The offline mode will analyze all images in a given directory that are 
taken by the LCOEngineering proposal (default), that are EXPOSE or EXPERIMENTAL
obstypes, and are defocused (FOCDMD != 0.0). To ensure that the proposal
is correctly included in the header, the Sequencer should be set to manual
when taking the images.

The real-time mode has the same criteria for which images to analyze, but
watches the input directory for new files. Right before an M2 tip/tilt run,
this watcher should be started on the directory where the files will appear. 
For the sbig cameras, this can be the `raw` directory for the given DAY-OBS.
The Sinistro frames need the amplifiers to be stitched together before they
can be analyzed, so the watcher should be directed to the `preview` directory
instead.

[//]: # (Code examples of the usage)
## Examples
The standard offline analysis can be run using the following command.
```
docker pull docker.lco.global/pupe-pat:{0.1.0}
docker run --rm -it -v /home/eng/pupe-pat-test-data:/input -v /home/eng/pupe-pat-tests-{0.1.0}:/output docker.lco.global/pupe-pat:{0.1.0} run_pupepat --output-dir /output --input-dir /input
```
Again replace `{0.1.0}` with the current version. This gets the most 
recent version and starts up a docker container
(which is kind of like a virtual machine).
The `--rm` flag removes the container when it is finished to keep the system
cleaned up. The `-it` runs the container in interactive mode. This should 
not be strictly necessary. The `-v` maps a directory on the current machine
into a directory in the docker container. The path on the left should be wherever
the raw data is and where you want the output results to go. The argument on the 
right of the `:` can be anything as long as the `--input-dir` and `output-dir`
are updated accordingly.

The only other parameters are `--output-table` and `--output-plot`. These control the
names of the output files produced at the end of the analysis. For example, if you 
wanted to run the analysis a second time but not overwrite the results, you could
run
```
docker pull docker.lco.global/pupe-pat:{0.1.0}
docker run --rm -it -v /home/eng/pupe-pat-test-data:/input -v /home/eng/pupe-pat-tests-{0.1.0}:/output docker.lco.global/pupe-pat:{0.1.0} run_pupepat --output-dir /output --input-dir /input --output-table pupe-pat-run2.dat --output-plot pupe-pat-run2
```
Note that the `--output-table` flag needs an extension (e.g. `.dat` or `.txt`),
but the `--output-plot` flag does not. It automatically makes pdf files. This difference
is because each focus position has a separate plot.

The real-time analysis can be run in the same way (all of the arguments are 
identical), but using a different command. 
```
docker pull docker.lco.global/pupe-pat:{0.1.0}
docker run --rm -it -v /home/eng/pupe-pat-test-data:/input -v /home/eng/pupe-pat-tests-{0.1.0}:/output docker.lco.global/pupe-pat:{0.1.0} run_pupepat_listener --output-dir /output --input-dir /input
```
Once images are done being taken for the night, press `control-c` to stop the listener
which will also generate the merged pdfs described above.


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
