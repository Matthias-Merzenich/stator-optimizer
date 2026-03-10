# Stator Optimizer
stator_optimizer.py is a program for optimizing the stator (unchanging part) of patterns in Life-like cellular automata. It is intended to be a faster, more memory-efficient, and easier to use replacement for Jeremy Dover's [CGOL Stator Reducer](https://github.com/jeremydover/CGOL-Stator-Reducer).

-Matthias Merzenich

## Installation
This program requires the following Python packages:
* `numpy`
* `ortools`
* `python-lifelib`
* `scipy`

You can install the dependencies using the provided `requirements.txt`:

```
pip install -r requirements.txt
```

The `lifelib` package requires a C++ compiler (gcc or clang). The initial run of stator_optimizer.py for a particular cellular automaton rule will take some time to compile code. Subsequent runs with the same rule will skip this step.

`lifelib` version 2.5.9 or higher is required when the input pattern is in rotor descriptor format. Since the version available on PyPI may be older, the program automatically checks your installation. If it detects the rotor descriptor format and an incompatible `lifelib` version, it attempts to clone the [latest version](https://gitlab.com/apgoucher/lifelib) to the current working directory using `git`.

## Acknowledgements
Thanks to Jeremy Dover for the original idea and implementation that inspired this program.

## Usage
```
usage: stator_optimizer.py [-h] [-a LEFT RIGHT TOP BOTTOM [{off,any}]]
                           [-d {rotor,stator} DISTANCE [{off,any}]]
                           [-s SYMMETRY] [-o [OBJECTIVES ...]]
                           [-c LEFT RIGHT TOP BOTTOM] [--solution_only]
                           input_file ticks

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
                           6 (gray)      | off         | either on or off
                        This only applies to stator cells. If a cell is
                        determined to be part of the rotor, then its state
                        will not be changed in the solution.
                        Accepted input formats are RLE (.rle), macrocell
                        (.mc), and rotor descriptor format.

  ticks                 Number of time steps to run the pattern for analysis.

options:
  -h, --help            Show this help message and exit.
  -a, --adjust LEFT RIGHT TOP BOTTOM [{off,any}]
                        Expand the search box by the given distances; negative
                        values contract it. An optional boundary state can
                        also be specified; if set to 'any', then the CA rules
                        will not be applied at the boundary.
  -d, --distance {rotor,stator} DISTANCE [{off,any}]
                        Force stator cells to be within the given distance of
                        either the rotor or the initial stator using the
                        taxicab metric. To set distances from both the rotor
                        and the stator, specify the option twice. For still
                        life searches, the distance is measured from the
                        forced cells, rather than the rotor. An optional
                        boundary state can also be specified; if set to 'any',
                        then the CA rules will not be applied at the boundary.
  -s, --symmetry SYMMETRY
                        Force the stator to have the given symmetry type.
                        Available types:
                          "C1" (default), "C2", "C4", "D2-", "D2|",
                          "D2/", "D2\", "D4+", "D4x", "D8".
                        Symmetry is applied relative to the center of the
                        search area. Always enclose the symmetry type in
                        double quotes (""); for D2| and D2\ symmetries you may
                        need to type "D2\|" and "D2\\", respectively.
  -o, --optimize [OBJECTIVES ...]
                        A list of properties to optimize ordered by priority.
                        Available properties:
                          min_pop, max_pop, area, width, height, symmetry,
                          left, right, top, bottom, nw, ne, sw, se, diag,
                          back_diag, min_change, max_change, boundary_pop,
                          valid_boundary.
                        The calculations of these properties only use the
                        stator. Append '_with_rotor' to also use the clipped
                        rotor from the box specified by the -c option.
  -c, --clip_rotor LEFT RIGHT TOP BOTTOM
                        Expand or contract the clipped rotor box. The initial
                        box is the bounding box of the input pattern. This
                        option is only used when '_with_rotor' objectives are
                        specified by the -o option.
  --solution_only       Only print the solution. If no solution exists or if
                        the input pattern is already optimal, print nothing.
```
