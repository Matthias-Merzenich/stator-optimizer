# Stator Optimizer
This is a simple Python script for reducing the stator (unchanging part) of
patterns in Life-like cellular automata. It is inspired by a similar program
by Jeremy Dover. At the moment, it has essentially no features and only works
in Conway's Game of Life.

This script requires the following Python packages:
* numpy
* python-lifelib
* ortools

-Matthias Merzenich

------------------------------------------------------------------------------
```
usage: stator_optimizer.py [-h] input_file ticks

positional arguments:
  input_file  Input LifeHistory file in RLE format
  ticks       Number of ticks to advance the pattern

options:
  -h, --help  show this help message and exit
```
