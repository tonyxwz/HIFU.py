import unittest
from pyHIFU.geometric import *

class Test_LinesMethods(unittest.TestCase):

    def test_upper(self):
        self.assertEqual('foo'.upper(), 'FOO')

    def test_isupper(self):
        self.assertTrue('FOO'.isupper())
        self.assertFalse('Foo'.isupper())

    def test_split(self):
        s = 'hello world'
        self.assertEqual(s.split(), ['hello', 'world'])
        # check that s.split fails when the separator is not a string
        with self.assertRaises(TypeError):
            s.split(2)


class Test_SurfacesMethods(unittest.TestCase):
    pass

class Test_VolumesMethods(unittest.TestCase):
    pass

if __name__ == '__main__':
    unittest.main()