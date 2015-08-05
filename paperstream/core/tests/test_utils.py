import numpy as np

from django.test import TestCase
from django.conf import settings
from django.core.exceptions import ValidationError

from ..utils import pad_vector, pad_neighbors

class PaddingTest(TestCase):

    def setUp(self):
        self.size = 10
        self.size2 = 15

    def test_padding_vector(self):
        vec = np.random.rand(self.size)
        vec_pad = pad_vector(vec)
        self.assertTrue((vec_pad[:self.size] == vec[:self.size]).all())
        self.assertTrue(all(x is None for x in vec_pad[self.size:]))

    def test_padding_vector_returns_list(self):
        vec = np.zeros((self.size, ))
        vec_pad = pad_vector(vec)
        self.assertTrue(isinstance(vec_pad, list))

    def test_padding_vector_can_take_nparray_as_input(self):
        vec = np.zeros((self.size, ))
        vec_pad = pad_vector(vec)
        self.assertTrue(isinstance(vec_pad, list))
        self.assertTrue((vec_pad[0:self.size] == vec[0:self.size]).all())
        self.assertTrue(all(x is None for x in vec_pad[self.size:]))

    def test_vector_canNOT_be_larger_than_maximum(self):
        vec = np.zeros((settings.NLP_MAX_VECTOR_SIZE+1, ))
        with self.assertRaises(ValidationError):
            pad_vector(vec)

    def test_vector_can_be_equal_to_maximum(self):
        vec = np.zeros((settings.NLP_MAX_VECTOR_SIZE, ))
        pad_vector(vec)

    def test_padding_neighbors(self):
        vec = np.random.rand(self.size)
        vec_pad = pad_neighbors(vec)
        self.assertTrue((vec_pad[0:self.size] == vec[0:self.size]).all())
        self.assertTrue(all(x is None for x in vec_pad[self.size:]))

    def test_padding_neighbors_returns_list(self):
        vec = np.zeros((self.size, ))
        vec_pad = pad_neighbors(vec)
        self.assertTrue(isinstance(vec_pad, list))

    def test_padding_neighbors_can_take_nparray_as_input(self):
        vec = np.zeros((self.size, ))
        vec_pad = pad_neighbors(vec)
        self.assertTrue(isinstance(vec_pad, list))
        self.assertTrue((vec_pad[0:self.size] == vec[0:self.size]).all())
        self.assertTrue(all(x is None for x in vec_pad[self.size:]))

    def test_vector_can_be_larger_than_maximum(self):
        vec = np.zeros((settings.NLP_MAX_KNN_NEIGHBORS + 1, ))
        with self.assertRaises(ValidationError):
            pad_neighbors(vec)
