import numpy as np
import tempfile
from unittest import skip

from test.test_pupil_synth import PupilSynthTestCase

from src.pupil_plate_image import PupilPlateImage
from src.source import Source, PointSource, CircularPupil


class TestPupilPlateImage(PupilSynthTestCase):
    def test___init__(self):
        ppi = PupilPlateImage()  # defaults to (4096, 4096)
        self.assertIsInstance(ppi, PupilPlateImage)
        image = ppi.render()
        self.assertEqual(image.shape, (4096, 4096))

    def test_add_source(self):
        ppi = PupilPlateImage((0, 0))
        ppi.add_source(PointSource(0, 0), 0, 0)
        self.assertEqual(len(ppi.sources), 1)

    def test_add_circular_pupil(self):
        ppi = PupilPlateImage()  # defaults to (4096, 4096)
        ppi.add_poisson_noise = False
        ppi.sky_level = 0
        ppi.read_noise = 0

        inner_radius, outer_radius, counts = 100, 200, 1000
        cp = CircularPupil(inner_radius, outer_radius, counts)

        # add the pupil to the center of the plate
        cx = cy = int(ppi.render().shape[0] / 2)
        ppi.add_source(cp, cx, cy)

        # check that the plate knows about it's new Source
        self.assertTrue(cp in ppi.sources)

        # check that the Source now knows where it is on the plate
        self.assertEqual((cx, cy), ppi.sources[0].position)

        # check that the Source is actually there.
        annulus_px = cp.x + cp.inner_radius + 10
        self.assertEqual(0, ppi.render()[cx, cy])  # center is 0
        self.assertEqual(cp.counts, ppi.render()[annulus_px, annulus_px])

    def test_add_point_source(self):
        ppi = PupilPlateImage()  # defaults to (4096, 4096)
        ppi.add_poisson_noise = False
        ppi.sky_level = 0
        ppi.read_noise = 0

        radius, counts = 10, 1000
        ps = PointSource(radius, counts)
        image = ppi.render()

        # add the pupil to the center of the plate
        cx = cy = int(image.shape[0] / 2)
        ppi.add_source(ps, cx, cy)

        # must re-render the image if a ppi parameter (like sources[]) has changed
        image = ppi.render()

        # check that the plate knows about it's new Source
        self.assertTrue(ps in ppi.sources)

        # check that the Source now knows where it is on the plate
        self.assertEqual((cx, cy), ppi.sources[0].position)

        # check that the Source is actually there.
        self.assertEqual(ps.counts, image[ps.x, ps.y])

    def test_pixel_perfect_placement(self):
        ppi = PupilPlateImage((4, 4))  # typical even number of pixels/edge
        ppi.add_poisson_noise = False
        ppi.sky_level = 0
        ppi.read_noise = 0

        radius, counts = 1, 99  # radius 0 gives single pixel decal
        ps = PointSource(radius, counts)

        # add the single pixel point source to the plate
        cx = 2
        cy = 1
        ppi.add_source(ps, cx, cy)

        image = ppi.render()

        # check that source is where we placed it
        for x in range(image.shape[0]):
            for y in range(image.shape[1]):
                if not ((x == cx) and (y == cy)):
                    self.assertEqual(0, image[x, y], '({x:}, {y:})'.format(x=x, y=y))
        self.assertEqual(counts, image[cx, cy], '({x:}, {y:})'.format(x=cx, y=cy))


class TestPupilPlateImageCutoutExtraction(PupilSynthTestCase):

    def setUp(self):
        self.ppi = PupilPlateImage(shape=(1024, 1024))
        self.cp = CircularPupil(100, 200, 1000.)
        self.ppi.add_source(self.cp, 512, 512)

    def test_cutout_margin_factor(self):
        expected_default = 1.2  # from the class definition
        # check the default value
        self.assertEqual(expected_default, self.ppi.cutout_margin_factor)

        # check setting/getting
        new_value = 1.618
        self.ppi.cutout_margin_factor = new_value
        self.assertEqual(new_value, self.ppi.cutout_margin_factor)

    def test_cutout_for_source_cutout(self):
        # if we set the cutout_margin_factor to 1.0,
        self.ppi.cutout_margin_factor = 1.0
        # then the shape of the decal and cutout should be equal
        cutout = self.ppi.cutout_for_source(self.cp)
        self.assertEqual(self.cp.decal.shape, cutout.shape)


#@skip
class TestPupilPlateImageFileIO(PupilSynthTestCase):
    """It might seem like we're just testing astropy.io.fits here,
    but I think we'll try to write the synthesized pupil parameters
    into the header pretty soon. That prepares us to run a FITS file
    through PUPE-PAT and know it's own correct answers.
    """

    def setUp(self):
        self.ppi = PupilPlateImage(shape=(1024, 1024))
        self.cp = CircularPupil(100, 200, 1000.)
        self.ppi.add_source(self.cp, 512, 512)

    def test_write_and_read(self):
        # create a temporary file with a name
        temp_file = tempfile.NamedTemporaryFile(dir='/tmp')
        temp_file.close()

        # write the PupilPlateImage to the file
        self.ppi.write_fits_to(temp_file.name)
        # read it into a new PPI instance
        ppi_from_file = PupilPlateImage()
        ppi_from_file.read_fits_from(temp_file.name)

        image_from_file = ppi_from_file.render()
        image_from_original = self.ppi.render()

        # compare the original and input PupilPlateImages
        self.assertEqual(image_from_original.all(), image_from_file.all())  # big

        # compare the header data with the Source parameters
        # TODO: implement FITS header data I/O
        pass


class TestPoissonMethod(PupilSynthTestCase):
    def test_poisson_methods(self):
        seed = 1234
        rs1 = np.random.RandomState(seed)
        rs2 = np.random.RandomState(seed)

        shape = (64, 64)
        sky = 100

        # method 1
        image1 = np.zeros(shape)
        image1 += sky
        image1 = rs1.poisson(image1)

        # method 2
        image2 = np.zeros(shape)
        shot_noise = rs2.poisson(sky, shape)
        image2 += shot_noise

        self.assertEqual(image1.all(), image2.all())  # looks like the methods are equivalent



