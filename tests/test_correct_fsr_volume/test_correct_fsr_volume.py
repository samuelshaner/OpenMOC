#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import PinCellInput


class CorrectFSRVolumeTestHarness(TestHarness):
    """Tests Geometry::correctFSRVolume() method."""

    def __init__(self):
        super(CorrectFSRVolumeTestHarness, self).__init__()
        self.input_set = PinCellInput()
        self._result = ''

    def _generate_tracks(self):
        """Generate Tracks and segments and correct the FSR volume for
        the fuel."""

        # Always use 1 thread for FSR reproducibility
        self.track_generator.setNumThreads(1)
        self.track_generator.generateTracks()

        old_volume = self.track_generator.getFSRVolume(1)
        self.track_generator.correctFSRVolume(1, 2.7)
        one_thread_volume = self.track_generator.getFSRVolume(1)

        # Save results for one thread
        self._result += str(round(old_volume,6)) + '\n'
        self._result += str(round(one_thread_volume, 1)) + '\n'

        # get results with multiple threads
        self.track_generator.setNumThreads(2)
        self.track_generator.correctFSRVolume(1, 2.5)
        two_thread_volume = self.track_generator.getFSRVolume(1)

        # Save results for two threads
        self._result += str(round(two_thread_volume,1))

    def _get_results(self, num_iters=False, keff=False, fluxes=False,
                     num_fsrs=True, num_segments=True, num_tracks=True,
                     hash_output=False):
        """Return the result string"""
        return self._result


if __name__ == '__main__':
    harness = CorrectFSRVolumeTestHarness()
    harness.main()
