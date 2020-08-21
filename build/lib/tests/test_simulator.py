#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import unittest
import numpy as np

sys.path.insert(0, '../')
import simulator


class TestSimulator(unittest.TestCase):
    def set_variables(self):
        self.mean_size = 1
        self.doubling_time = 18
        self.div_steps = 10
        self.ncells = 5000
        self.gr = np.log(2) / self.doubling_time
        self.sim = simulator.Simulator(ncells=self.ncells, gr=self.gr, sb=self.mean_size, steps=self.div_steps)

    def test_newgr(self):
        self.set_variables()
        self.assertEqual(self.sim.newgr(0), 1)
        for i in range(1, 6):
            self.assertNotEqual(self.sim.newgr(i), 1)
        for i in np.arange(0.1, 1, 0.1):
            self.assertNotEqual(self.sim.newgr(i), 1)
        for i in np.arange(0.01, 1, 0.01):
            self.assertNotEqual(self.sim.newgr(i), 1)

    def test_newdivpar(self):
        self.set_variables()
        self.assertEqual(self.sim.newdivpar(0), 0.5)
        for i in np.arange(0.1, 1, 0.1):
            self.assertNotEqual(self.sim.newdivpar(i), 0.5)
        for i in np.arange(0.01, 0.9, 0.01):
            self.assertNotEqual(self.sim.newdivpar(i), 0.5)

    

    def test_get_sz(self):
        self.set_variables()
        for i in range(len(self.sim.cells)):
            self.assertEqual(self.sim.get_sz(i), self.sim.cells[i].V)

    def test_get_ndiv(self):
        self.set_variables()
        for i in range(len(self.sim.cells)):
            self.assertEqual(self.sim.get_ndiv(i), self.sim.cells[i].ndiv)

    def test_get_gr(self):
        self.set_variables()
        for i in range(len(self.sim.cells)):
            self.assertEqual(self.sim.get_gr(i), self.sim.cells[i].gr)

    def test_get_dp(self):
        self.set_variables()
        for i in range(len(self.sim.cells)):
            self.assertEqual(self.sim.get_dp(i), self.sim.cells[i].dp)

    def test_get_next_t(self):
        self.set_variables()
        for i in range(len(self.sim.cells)):
            self.assertEqual(self.sim.get_next_t(i), self.sim.cells[i].nextt)


if __name__ == "__main__":
    unittest.main()