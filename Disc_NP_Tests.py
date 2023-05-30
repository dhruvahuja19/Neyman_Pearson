import unittest
import Neyman_Pearson as NP

class MyTestCase(unittest.TestCase):
    def test_F20_Q9(self):
        nullhypothesis = [(i, 0.05) for i in range(1, 21)]
        alternative = [(i, 0.025) for i in range(1, 11)] + [(i, 0.075) for i in range(11, 21)]
        alpha = 0.05
        vals = NP.disc_NP(nullhypothesis, alternative, alpha)
        print(vals[0][1])
        self.assertEqual(vals[0][1], [[1, 2, 3, 4, 5, 6, 7,  8,  9, 10]])
        self.assertEqual(vals[1][1], [])
        self.assertEqual(vals[2][1], [[11, 12, 13, 14, 15, 16, 17, 18, 19, 20]])
    def test_sp20_Q9(self):
        nullhypothesis = [(i, 0.25) for i in range(1, 5)]
        alternative = [(1, 0.25), (2, (1/6)),(3, (1.0/2)), (4, (1.0/12))]
        alpha = (0.3)
        vals = NP.disc_NP(nullhypothesis, alternative, alpha)






if __name__ == '__main__':
    unittest.main()
