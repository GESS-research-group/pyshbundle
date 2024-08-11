import unittest
from .validation_pyshbundle import validation_pyshbundle

class TestpyshbundleAccuracy(unittest.TestCase):
    
    def test_pyshbundle(self):
        result = validation_pyshbundle() 
        expected = "expected_result"
        self.assertEqual(result, expected)

if __name__ == '__main__':
    unittest.main()