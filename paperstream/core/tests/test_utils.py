import numpy as np

from django.test import TestCase
from django.conf import settings
from django.core.exceptions import ValidationError

from ..utils import pad_vector, pad_neighbors

class PaddingTest(TestCase):

    def setUp(self):
        settings.NLP_MAX_VECTOR_SIZE = 7
        settings.NLP_MAX_KNN_NEIGHBORS = 6

    def test_padding_vector(self):
        vec = [1., 2., 3., 4., 5.]
        vec_pad = pad_vector(vec)
        self.assertTrue(vec_pad[0:5] == vec[0:5])
        self.assertTrue(vec_pad[5:] == [None, None])

    def test_padding_vector_returns_list(self):
        vec = np.zeros((5, ))
        vec_pad = pad_vector(vec)
        self.assertTrue(isinstance(vec_pad, list))

    def test_padding_vector_can_take_nparray_as_input(self):
        vec = np.zeros((5, ))
        vec_pad = pad_vector(vec)
        self.assertTrue(isinstance(vec_pad, list))
        self.assertTrue((vec_pad[0:5] == vec[0:5]).all())
        self.assertTrue(vec_pad[5:] == [None, None])

    def test_vector_can_be_larger_than_maximum(self):
        vec = np.zeros((settings.NLP_MAX_VECTOR_SIZE+1, ))
        with self.assertRaises(ValidationError):
            pad_vector(vec)

    def test_padding_neighbors(self):
        vec = [1., 2., 3., 4., 5.]
        vec_pad = pad_neighbors(vec)
        self.assertTrue(vec_pad[0:5] == vec[0:5])
        self.assertTrue(vec_pad[5:] == [None])

    def test_padding_neighbors_returns_list(self):
        vec = np.zeros((5, ))
        vec_pad = pad_neighbors(vec)
        self.assertTrue(isinstance(vec_pad, list))

    def test_padding_neighbors_can_take_nparray_as_input(self):
        vec = np.zeros((5, ))
        vec_pad = pad_neighbors(vec)
        self.assertTrue(isinstance(vec_pad, list))
        self.assertTrue((vec_pad[0:5] == vec[0:5]).all())
        self.assertTrue(vec_pad[5:] == [None])

    def test_vector_can_be_larger_than_maximum(self):
        vec = np.zeros((settings.NLP_MAX_KNN_NEIGHBORS + 1, ))
        with self.assertRaises(ValidationError):
            pad_neighbors(vec)
