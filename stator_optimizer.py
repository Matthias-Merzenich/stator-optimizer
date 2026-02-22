#!/usr/bin/env python3

import argparse
import re
import os

import lifelib
import numpy as np
from ortools.sat.python import cp_model

# Treat coordinate tuples as single objects
def coordinate_view(arr):
        return np.ascontiguousarray(arr).view(np.dtype((np.void, arr.dtype.itemsize * arr.shape[1])))

# For every cell in the numpy array origin_cells, count the number
# of its neighbours that are in the numpy array neighbour_cells.
def count_neighbours(origin_cells, neighbour_cells):
    neighbour_offsets = np.array([  [-1, -1], [-1, 0], [-1, 1],
                                    [ 0, -1],          [ 0, 1],
                                    [ 1, -1], [ 1, 0], [ 1, 1]
                                 ])

    # New array containing all possible neighbours of the
    # origin_cells arranged in groups of 8 by origin cell.
    all_neighbours = origin_cells[:, np.newaxis, :] + neighbour_offsets
    all_neighbours = all_neighbours.reshape(-1, 2)

    # Identify which cells are in neighbour_cells and sum them.
    is_neighbour = np.isin(coordinate_view(all_neighbours), coordinate_view(neighbour_cells))
    neighbour_counts = is_neighbour.reshape(-1, 8).sum(axis=1)

    return neighbour_counts

def get_rule(rle):
    rulestring = lifelib.load_rules("b3s23").lifetree().pattern(rle).getrule()

    if rulestring.endswith("History"):
        rulestring = rulestring[:-len("History")]

    match = re.match(r"^b([2-8]*)s([0-8]*)$", rulestring)
    assert match, "Invalid rule: only Life-like rules without B0 or B1 are supported"
    birth_digits, survival_digits = match.groups()
    assert survival_digits, "Invalid rule: no survival conditions"

    birth_at = [False] * 9
    survival_at = [False] * 9
    for digit in birth_digits:
        birth_at[int(digit)] = True
    for digit in survival_digits:
        survival_at[int(digit)] = True

    return rulestring, birth_at, survival_at

def clean_rle(rle):
    return "\n".join([line for line in rle.splitlines() if not line.startswith('#')])

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawTextHelpFormatter,
    description="A program to optimize the stator of patterns in Life-like cellular automata.",
    add_help=False
)
parser.add_argument(
    'input_file',
    help=   "File containing the pattern to be optimized. When the\n"
            "pattern is saved using a history rule, you can specify\n"
            "stator cell states in the solution as follows:\n"
            "------------------------------------------------\n"
            " history state | input state | solution state   \n"
            "---------------+-------------+------------------\n"
            " 0 (black)     | off         | either on or off \n"
            " 1 (green)     | on          | either on or off \n"
            " 2 (blue)      | off         | off              \n"
            " 3 (white)     | on          | off              \n"
            " 4 (red)       | off         | on               \n"
            " 5 (yellow)    | on          | on               \n"
            "This only applies to stator cells. If a cell is\n"
            "determined to be part of the rotor, then its state\n"
            "will not be changed in the solution.\n\n"
)
parser.add_argument(
    'ticks',
    type=int,
    help="Number of time steps to run the pattern for analysis."
)
parser.add_argument(
    '-h', '--help',
    action='help',
    default=argparse.SUPPRESS,
    help='Show this help message and exit.'
)
parser.add_argument(
    '-a', '--adjust',
    type=int,
    nargs=4,
    default=[0, 0, 0, 0],
    metavar=("LEFT", "RIGHT", "TOP", "BOTTOM"),
    help=   "Expand the search box by the given distances. Negative\n"
            "values contract the search box."
)
args = parser.parse_args()

ticks = args.ticks
adjust_left, adjust_right, adjust_top, adjust_bottom = args.adjust

if os.path.isfile(args.input_file):
    with open(args.input_file, 'r') as file:
        rle = file.read()
else:
    assert False, "Pattern file not found"

rulestring, birth_at, survival_at = get_rule(rle)

sess = lifelib.load_rules(rulestring, "bs8", "b12345678s012345678")
lt = sess.lifetree(n_layers=1)
lt4 = sess.lifetree(n_layers=4)

# Convert the input file to a history RLE. It's faster to change the
# cell states by modifying the RLE than to change them in the lifetree.
initial_pattern = lt4.pattern("", rulestring + "History")
initial_pattern += lt4.pattern(rle)
rle = initial_pattern.rle_string()

# Identify cells to be forced on or off in the solution
rle_forced_on = rle.translate(str.maketrans({   'A': 'B',
                                                'C': 'B',
                                                'D': 'E',
                                                'F': 'B'
                                            }))
rle_forced_off = rle.translate(str.maketrans({  'A': 'B',
                                                'B': 'E',
                                                'C': 'E',
                                                'D': 'B',
                                                'E': 'B',
                                                'F': 'B'
                                             }))
forced_on = lt4.pattern(rle_forced_on).layers()[2]
forced_off = lt4.pattern(rle_forced_off).layers()[2]

# I was getting segmentation faults when forced_on or
# forced_off were empty, so I had to include these checks.
# I don't know why this occurs.
if forced_on.nonempty():
    forced_on = lt.pattern("", rulestring) + forced_on
else:
    forced_on = lt.pattern("", rulestring)
if forced_off.nonempty():
    forced_off = lt.pattern("", rulestring) + forced_off
else:
    forced_off = lt.pattern("", rulestring)

# Convert all on cells to state 5 and all off cells to state 2.
# In LifeHistory, State 5 cells will only remain state 5 as long as
# they never die, so they can be used to detect the initial stator.
rle = rle.translate(str.maketrans({ 'A': 'E',
                                    'C': 'E',
                                    'D': 'B',
                                    'F': 'B'
                                  }))

initial_pattern = lt4.pattern(rle)
assert initial_pattern.nonempty(), "Input pattern is empty"
_, _, width, height = initial_pattern.bounding_box

# Remove state 2 cells
mask = lt4.pattern("", rulestring + "History")
mask[0:width, 0:height] = 2
initial_pattern -= mask

initial_population = initial_pattern.population
print("\nAnalyzing pattern...")
print(f'Initial population: {initial_population}')
print(f'Initial bounding box: {width} x {height}')
print(f'Search bounding box: {width + adjust_left + adjust_right} x {height + adjust_top + adjust_bottom}')

final_pattern = initial_pattern[ticks]
envelope = final_pattern.layers()[1]
initial_stator_on = final_pattern.layers()[2]
rotor = lt.pattern("", "bs8")
rotor += envelope - initial_stator_on

# Restrict the stator to the search box
# with an additional 1-cell-thick border.
stator = lt.pattern("", "b12345678s012345678")
stator[     0 - adjust_left : width + adjust_right,
            0 - adjust_top : height + adjust_bottom
      ] = 1
stator -= rotor
forced_on &= stator
forced_off &= stator
boundary_stator = stator[1] - stator - rotor
stator += boundary_stator

# Rotor cells that have a stator cell as a neighbour. Note that stator[1]
# comes first in the assignment, because we want adjacent_rotor to have
# rule B12345678/S012345678 for later manipulation.
adjacent_rotor = stator[1] & (rotor - rotor[1])

initial_pattern_2_state = lt.pattern("", rulestring)
initial_pattern_2_state += initial_pattern.layers()[2]
rotor_phases = [adjacent_rotor[1] & rotor & initial_pattern_2_state[t] for t in range(ticks + 1)]

# For each generation calculate the envelope covering all cells
# that changed in the last tick or that border a changed cell.
# change_envelopes[0] assumes all adjacent rotor cells changed.
change_envelopes = [adjacent_rotor[1]]
for t in range(ticks):
    changed_cells = rotor_phases[t] ^ rotor_phases[t+1]
    change_envelopes.append(changed_cells[1] & adjacent_rotor[1])

# Coordinate lists of various groups of cells
stator_cells = stator.coords().tolist()
stator_cells = set(map(tuple, stator_cells))
forced_on_cells = forced_on.coords().tolist()
forced_on_cells = set(map(tuple, forced_on_cells))
forced_off_cells = (forced_off + boundary_stator).coords().tolist()
forced_off_cells = set(map(tuple, forced_off_cells))
non_adjacent_stator_cells = (stator - adjacent_rotor[1]).coords().tolist()

# Set of tuples (x,y,c) where (x,y) is a stator cell and c is a count of
# its on rotor neighbours in some generation. A single stator cell will
# correspond to multiple elements in this set.
stator_neighbour_counts = set()
for t in range(ticks):
    gen_t_rotor_coords = rotor_phases[t].coords()
    gen_t_stator_coords = (change_envelopes[t] & stator).coords()

    gen_t_neighbour_counts = count_neighbours(gen_t_stator_coords, gen_t_rotor_coords)
    gen_t_neighbour_counts = np.expand_dims(gen_t_neighbour_counts, axis=1)

    coords_and_counts = np.hstack(( gen_t_stator_coords,
                                    gen_t_neighbour_counts
                                  )).astype(np.int64).tolist()
    stator_neighbour_counts.update(set(map(tuple, coords_and_counts)))

for x,y in non_adjacent_stator_cells:
    stator_neighbour_counts.add((x,y,0))

# Set of tuples (x, y, c, state_0, state_1) where (x,y) is a rotor cell, c is a
# count of its on rotor neighbours in some generation, state_0 is the state of
# the cell in that generation, and state_1 is the state in the next generation.
rotor_transitions = set()
for t in range(ticks):
    gen_t_rotor_coords = rotor_phases[t].coords()
    gen_t_adjacent_rotor_coords = (change_envelopes[t] & adjacent_rotor).coords()

    gen_t_rotor_state_0 = rotor_phases[t][gen_t_adjacent_rotor_coords]
    gen_t_rotor_state_1 = rotor_phases[t+1][gen_t_adjacent_rotor_coords]
    gen_t_rotor_neighbour_counts = count_neighbours(gen_t_adjacent_rotor_coords, gen_t_rotor_coords)

    gen_t_rotor_state_0 = np.expand_dims(gen_t_rotor_state_0, axis=1)
    gen_t_rotor_state_1 = np.expand_dims(gen_t_rotor_state_1, axis=1)
    gen_t_rotor_neighbour_counts = np.expand_dims(gen_t_rotor_neighbour_counts, axis=1)

    coords_and_counts = np.hstack(( gen_t_adjacent_rotor_coords,
                                    gen_t_rotor_neighbour_counts,
                                    gen_t_rotor_state_0,
                                    gen_t_rotor_state_1
                                  )).astype(np.int64).tolist()
    rotor_transitions.update(set(map(tuple, coords_and_counts)))

print("\nSetting up optimization search...")

model = cp_model.CpModel()

stator_int = {}
stator_bool = {}
for x,y in stator_cells:
    stator_int[x,y] = model.NewIntVar(0,1,f'stator_int_{x}_{y}')
    stator_bool[x,y] = model.NewBoolVar(f'stator_bool_{x}_{y}')
    model.Add(stator_int[x,y] == 1).OnlyEnforceIf(stator_bool[x,y])
    model.Add(stator_int[x,y] == 0).OnlyEnforceIf(stator_bool[x,y].Not())

    if (x,y) in forced_off_cells:
        model.Add(stator_int[x,y] == 0)
    elif (x,y) in forced_on_cells:
        model.Add(stator_int[x,y] == 1)

# Apply the CA rules to the stator cells
for x,y,rotor_sum in stator_neighbour_counts:
    stator_sum = sum(   stator_int[x+u, y+v]
                        for u in [-1,0,1] for v in [-1,0,1]
                        if (u,v) != (0,0) and (x+u, y+v) in stator_cells
                    )
    neighbour_sum = rotor_sum + stator_sum
    if rulestring == "b3s23":
        model.AddLinearConstraint(neighbour_sum, 2, 3).OnlyEnforceIf(stator_bool[x,y])
        model.Add(neighbour_sum != 3).OnlyEnforceIf(stator_bool[x,y].Not())
    else:
        for count in range(9):
            if birth_at[count]:
                model.Add(neighbour_sum != count).OnlyEnforceIf(stator_bool[x,y].Not())
            if not survival_at[count]:
                model.Add(neighbour_sum != count).OnlyEnforceIf(stator_bool[x,y])

# Apply the CA rules to the rotor transitions
for x,y,rotor_sum,state_0,state_1 in rotor_transitions:
    stator_sum = sum(   stator_int[x+u, y+v]
                        for u in [-1,0,1] for v in [-1,0,1]
                        if (u,v) != (0,0) and (x+u, y+v) in stator_cells
                    )
    neighbour_sum = rotor_sum + stator_sum
    if rulestring == "b3s23":
        if state_0 == 0 and state_1 == 0:
            model.Add(neighbour_sum != 3)
        if state_0 == 0 and state_1 == 1:
            model.Add(neighbour_sum == 3)
        if state_0 == 1 and state_1 == 0:
            model.Add(neighbour_sum != 2)
            model.Add(neighbour_sum != 3)
        if state_0 == 1 and state_1 == 1:
            model.AddLinearConstraint(neighbour_sum, 2, 3)
    else:
        for count in range(9):
            if state_0 == 0 and state_1 == 0 and birth_at[count]:
                model.Add(neighbour_sum != count)
            if state_0 == 0 and state_1 == 1 and not birth_at[count]:
                model.Add(neighbour_sum != count)
            if state_0 == 1 and state_1 == 0 and survival_at[count]:
                model.Add(neighbour_sum != count)
            if state_0 == 1 and state_1 == 1 and not survival_at[count]:
                model.Add(neighbour_sum != count)

# Use the input pattern as a hint to the solver
intial_stator_on_cells = set(map(tuple, (stator & initial_stator_on).coords().tolist()))
intial_stator_off_cells = set(map(tuple, (stator - initial_stator_on).coords().tolist()))
for x,y in intial_stator_on_cells:
    model.AddHint(stator_int[x,y], 1)
for x,y in intial_stator_off_cells:
    model.AddHint(stator_int[x,y], 1)

size = sum(stator_int[x,y] for x,y in stator_cells)
model.Minimize(size)

solver = cp_model.CpSolver()

print("Beginning optimization search...")

status = solver.Solve(model)

if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
    if status == cp_model.OPTIMAL:
        print("\nOptimal solution found.")
    elif status == cp_model.FEASIBLE:
        print("\nNon-optimal solution found.")

    output_pattern = initial_pattern_2_state & rotor
    for x,y in stator_cells:
        if solver.Value(stator_int[x,y]) == 1:
            output_pattern[x,y] = 1

    # Only print the pattern if:
    # * it has a smaller population than the original, or
    # * the original stator contains a cell outside the search area, or
    # * there is a cell in the original stator that is forced off in the solution, or
    # * there is a cell not in the original stator that is forced on in the solution.
    if (
        output_pattern.population < initial_population or
        (initial_stator_on - (stator - boundary_stator)).nonempty() or
        (forced_off & initial_stator_on).nonempty() or
        (forced_on - initial_stator_on).nonempty()
    ):
        print(f"\nFinal population: {output_pattern.population}")
        print(clean_rle(output_pattern.rle_string()))
    elif status == cp_model.OPTIMAL:
        print("\nInput pattern is already optimal.")
    else:
        print("\nNew pattern is not smaller than the input pattern.")
elif status == cp_model.UNKNOWN:
    print("\nUNKNOWN: no solutions were found, but the model was not proven infeasible.")
elif status == cp_model.INFEASIBLE:
    print("\nINFEASIBLE: no solution possible.")

print()
