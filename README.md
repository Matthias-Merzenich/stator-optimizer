# Stator Optimizer
This is a simple Python script for optimizing the stator (unchanging part) of patterns in Life-like cellular automata. It is intended to be a faster, more memory-efficient replacement for Jeremy Dover's [CGOL Stator Reducer](https://github.com/jeremydover/CGOL-Stator-Reducer).

This script requires the following Python packages:
* numpy
* ortools
* python-lifelib
* scipy

-Matthias Merzenich

------------------------------------------------------------------------------
```
usage: stator_optimizer.py [-h] [-a LEFT RIGHT TOP BOTTOM] input_file ticks

A program to optimize the stator of patterns in Life-like cellular automata.

positional arguments:
  input_file            File containing the pattern to be optimized. When the
                        pattern is saved using a history rule, you can specify
                        stator cell states in the solution as follows:
                        ------------------------------------------------
                         history state | input state | solution state
                        ---------------+-------------+------------------
                         0 (black)     | off         | either on or off
                         1 (green)     | on          | either on or off
                         2 (blue)      | off         | off
                         3 (white)     | on          | off
                         4 (red)       | off         | on
                         5 (yellow)    | on          | on
                        This only applies to stator cells. If a cell is
                        determined to be part of the rotor, then its state
                        will not be changed in the solution.

  ticks                 Number of time steps to run the pattern for analysis.

options:
  -h, --help            Show this help message and exit.
  -a LEFT RIGHT TOP BOTTOM, --adjust LEFT RIGHT TOP BOTTOM
                        Expand the search box by the given distances. Negative
                        values contract the search box.
```
