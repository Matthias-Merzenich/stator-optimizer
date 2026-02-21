# Stator Optimizer
This is a simple Python script for reducing the stator (unchanging part) of patterns in Life-like cellular automata. This project was inspired by Jeremy Dover's [CGOL Stator Reducer](https://github.com/jeremydover/CGOL-Stator-Reducer). It features significantly improved speed and memory use, but currently has fewer options.

This script requires the following Python packages:
* numpy
* python-lifelib
* ortools

-Matthias Merzenich

------------------------------------------------------------------------------

```
usage: stator_optimizer.py [-h] [-a LEFT RIGHT TOP BOTTOM] input_file ticks

positional arguments:
  input_file            File containing the pattern to be optimized
  ticks                 Number of time steps to advance the pattern for analysis

options:
  -h, --help            show this help message and exit
  -a LEFT RIGHT TOP BOTTOM, --adjust LEFT RIGHT TOP BOTTOM
                        Expand the search box by the given distances (negative values contract the search box)
```
