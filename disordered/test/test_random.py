import random
import unittest

class TestSequenceFunctions(unittest.TestCase):
		# def setUp(self):
		# 	self.seq = range(10)
		# 
		# def test_shuffle(self):
		# 	random.shuffle(self.seq)
		# 	self.seq.sort()
		# 	self.assertEqual(self.seq, range(10))
		# 	
		# 	# Should raise an exception for an immutable sequence
		# 	self.assertRaises(TypeError, random.shuffle, (1,2,3))
		# 
		# def test_choice(self):
		# 	element = random.choice(self.seq)
		# 	self.assertTrue(element in self.seq)
		# 
		# def test_sample(self):
		# 	with self.assertRaises(ValueError):
		# 		random.sample(self.seq, 20)
		# 		
		# 	for element in random.sample(self.seq, 5):
		# 		self.assertTrue(element in self.seq)
		
		# def runTest(self):
		# 	widget = Widget('The widget')
		# 	self.assertEqual(widget.size(), (50, 50), 'incorrect default size')

# if __name__ == '__main__':
# 	unittest.main()

suite = unittest.TestLoader().loadTestsFromTestCase(TestSequenceFunctions)
unittest.TextTestRunner(verbosity=2).run(suite)