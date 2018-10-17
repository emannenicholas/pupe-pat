import numpy as np
from pupepat.test.pupil-synth.test.test_pupil_synth import PupilSynthTestCase, getsize

from pupepat.test.src.source import Source, PointSource, CircularPupil, EllipticalPupil


class TestSource(PupilSynthTestCase):
    def test___init__(self):
        # raise TypeError upon instanciation of abc.ABC
        self.assertRaises(TypeError, Source, 0, 0)


class TestPointSource(PupilSynthTestCase):
    def test___init__(self):
        x, y, radius, counts = 0, 0, 10, 100
        ps = PointSource(x, y, radius, counts)

        self.assertIsInstance(ps, PointSource)
        print(ps)

    def test_position(self):
        x, y, radius, counts = 0, 0, 10, 100
        ps0 = PointSource(radius, counts, x, y)
        self.assertEqual((x, y), ps0.position)

        ps1 = PointSource(radius, counts)
        self.assertEqual((None, None), ps1.position)

    def test_decal_shape(self):
        radius, counts = 10, 100
        ps = PointSource(radius, counts)
        expected_shape = (2*radius+1, 2*radius+1)
        self.assertEqual(expected_shape, ps.decal.shape)
        # also test source.half_with is correct
        self.assertEqual(radius, ps.half_width)

    def test_decal_sum(self):
        radius, counts = 1, 100
        ps = PointSource(radius, counts)
        self.assertEqual(counts, ps.decal.sum())

        radius = 2
        ps = PointSource(radius, counts)
        self.assertEqual(9 * counts, ps.decal.sum())

    def test_memory_size_is_constant(self):
        ps100 = PointSource(100, 11)
        ps1000 = PointSource(1000, 33)
        # size of decal should not change size of object
        self.assertEqual(getsize(ps100), getsize(ps1000))


class TestCircularPupil(PupilSynthTestCase):
    def test___init__(self):
        x, y, inner_radius, outer_radius, counts = 0, 0, 100, 200, 100
        cp = CircularPupil(inner_radius, outer_radius, counts, x=x, y=y)

        self.assertIsInstance(cp, CircularPupil)
        print(cp)

    def test_decal_shape(self):
        """the decal is shaped on the outer_radius
        """
        x, y, inner_radius, outer_radius, counts = 0, 0, 100, 200, 100
        cp = CircularPupil(inner_radius, outer_radius, counts, x=x, y=y)

        expected_shape = (2*outer_radius+1, 2*outer_radius+1)
        self.assertEqual(expected_shape, cp.decal.shape)
        # also test source.half_with is correct
        self.assertEqual(outer_radius, cp.half_width)

    def test_corner_values(self):
        """We're testing withing the decal, so x, y (the data coords) don't matter.
        """
        x, y, inner_radius, outer_radius, counts = 0, 0, 100, 200, 100
        cp = CircularPupil(inner_radius, outer_radius, counts, x=x, y=y)

        corner = 2*outer_radius  # corner of the decal
        self.assertEqual(0, cp.decal[0, 0])
        self.assertEqual(0, cp.decal[0, corner])
        self.assertEqual(0, cp.decal[corner, corner])
        self.assertEqual(0, cp.decal[corner, 0])

    def test_center_value(self):
        x, y, inner_radius, outer_radius, counts = 0, 0, 100, 200, 100
        cp = CircularPupil(inner_radius, outer_radius, counts, x=x, y=y)

        ci = cj = outer_radius + 1  # the coords of the middle pixel
        self.assertEqual(0, cp.decal[ci, cj])

    def test_decal_sum(self):
        x, y, inner_radius, outer_radius, counts = 0, 0, 1, 2, 100
        cp = CircularPupil(inner_radius, outer_radius, counts, x=x, y=y)

        expected_sum = (9-1)*counts  # 9 from outer - 1 from inner
        self.assertEqual(expected_sum, cp.decal.sum())

    def test_raises_ValueError(self):
        """If the inner_radius and center_vector are not compatible with outer_radius, raise ValueError

        The disk of the inner radius must be contained in the outer_radius.
        """
        inner_radius, outer_radius, counts = 3, 2, 100
        self.assertRaises(ValueError, CircularPupil, inner_radius, outer_radius, counts)

    def test_memory_size_is_constant(self):
        ps100 = CircularPupil(100, 200, 1000)
        ps100_size = getsize(ps100)
        ps100a = CircularPupil(100, 200, 1000)
        ps100a_size = getsize(ps100a)
        ps1000 = CircularPupil(1000, 2000, 1)
        ps1000_size = getsize(ps1000)
        # magnitude of counts should not matter
        self.assertEqual(ps100_size, ps100a_size)

        # size of the decal should not matter
        self.assertEqual(ps100_size, ps1000_size)

        # size must not change if decal is generated
        decal = ps1000.decal
        self.assertEqual(ps1000_size, getsize(ps1000))


class TestEllipticalPupil(PupilSynthTestCase):
    def setUp(self):
        self.outer_a = 400
        self.outer_b = 360
        self.outer_phi = np.pi/3.
        self.inner_a = 200
        self.inner_b = 180
        self.inner_phi = np.pi/3.
        self.counts = 100

        self.ep0 = EllipticalPupil(
            self.inner_a, self.inner_b, 0.,
            self.outer_a, self.outer_b, 0.,
            self.counts)

        self.ep1 = EllipticalPupil(
            self.inner_a, self.inner_b, self.inner_phi,
            self.outer_a, self.outer_b, self.outer_phi,
            self.counts)

        self.ep2 = EllipticalPupil(
            self.inner_a, self.inner_b, self.inner_phi,
            self.outer_a, self.outer_b, self.outer_phi,
            self.counts, center_vector=(11, 14))

    def test_decal_shape(self):
        expected_shape = (2*self.outer_a + 1, 2*self.outer_a + 1)
        self.assertEqual(expected_shape, self.ep1.decal.shape)
        # also test source.half_with is correct
        expected_half_width = max(self.outer_a, self.outer_b)
        self.assertEqual(expected_half_width, self.ep2.half_width)

    def test_some_pixel_values(self):
        # the center of the decal should be zero (or inner ring if center_vector is non-zero)
        cx, cy = self.ep1.half_width, self.ep1.half_width
        self.assertEqual(0, self.ep1.decal[cx, cy])

        # the center of inner ring should be zero for non-zero center_vector
        cx = self.ep2.half_width + self.ep2.center_vector[0]
        cy = self.ep2.half_width + self.ep2.center_vector[1]
        self.assertEqual(0, self.ep1.decal[cx, cy])

        # the vertices and co-vertices of the outer ring should have counts
        cx, cy = self.ep0.half_width, self.ep0.half_width
        self.assertEqual(self.counts, self.ep0.decal[cx + self.outer_b - 1, cy])
        self.assertEqual(self.counts, self.ep0.decal[cx - self.outer_b + 1, cy])
        self.assertEqual(self.counts, self.ep0.decal[cx, cy + self.outer_a - 1])
        self.assertEqual(self.counts, self.ep0.decal[cx, cy - self.outer_a + 1])

    def test_a_less_than_b(self):
        # make a new EllipticalPupil with As and Bs swapped and compare to self.ep1
        x_outer_a = self.outer_b
        x_outer_b = self.outer_a
        x_inner_a = self.inner_b
        x_inner_b = self.inner_a
        x_ep = EllipticalPupil(x_inner_a, x_inner_b, self.inner_phi,
                               x_outer_a, x_outer_b, self.outer_phi, self.counts)

        print(self.ep1)
        print(x_ep)

        # decal shapes should be the same
        self.assertEqual(self.ep1.decal.shape, x_ep.decal.shape)

        # the sum of the pixels in each decal should be the same
        self.assertEqual(np.sum(self.ep1.decal), np.sum(x_ep.decal))

        # if we subtract the decals, the np.sum of the difference should be 0
        # b/c while most pixels become zero, the non-zero sums cancel
        self.assertEqual(0, np.sum(np.subtract(self.ep1.decal, x_ep.decal)))

        # also test source.half_with is correct
        expected_half_width = max(self.outer_a, self.outer_b)
        self.assertEqual(expected_half_width, self.ep2.half_width)

    def test_memory_size_is_constant(self):
        ia, ib, i_phi, oa, ob, o_phi = 400, 360, np.pi/5., 200, 180, np.pi/6
        counts = 11
        ep1 = EllipticalPupil(ia, ib, i_phi, oa, ob, o_phi, counts)
        ep2 = EllipticalPupil(ia-10, ib-10, i_phi, oa+10, ob, o_phi, counts)
        ep1_size = getsize(ep1)
        ep2_size = getsize(ep2)
        # different decal sizes should not matter
        self.assertEqual(ep1_size, ep2_size)

        # make sure that size stays the same after making the decal
        d1 = ep1.decal  # create the decal
        self.assertEqual(ep1_size, getsize(ep1))
