#!/usr/bin/env python3

import os
import argparse
import numpy as np
import lifelib
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

parser = argparse.ArgumentParser()
parser.add_argument(
    'input_file',
    help="Input LifeHistory file in RLE format"
)
parser.add_argument(
    'ticks',
    type=int,
    help="Number of ticks to advance the pattern"
)
args = parser.parse_args()

ticks = args.ticks

if os.path.isfile(args.input_file):
    with open(args.input_file, 'r') as file:
        rle = file.read()
else:
    assert False, "Search pattern file not found"

sess = lifelib.load_rules("b3s23", "bs8", "b12345678s012345678")
lt = sess.lifetree(n_layers=1)
lt4 = sess.lifetree(n_layers=4)

# Convert all on cells to state 5 and all off cells to state 2.
# In LifeHistory, State 5 cells will only remain state 5 as long as
# they never die, so they can be used to detect the initial stator.
rle = rle.translate(str.maketrans({ 'A': 'E',
                                    'C': 'E',
                                    'D': 'B',
                                    'F': 'B'
                                  }))

initial_pattern = lt4.pattern(rle)
_, _, width, height = initial_pattern.bounding_box

# Remove state 2 cells
mask = lt4.pattern("","LifeHistory")
mask[0:width, 0:height] = 2
initial_pattern -= mask

initial_population = initial_pattern.population
print("\nAnalyzing pattern...")
print(f'Initial population: {initial_population}')
print(f'Bounding box: {width} x {height}')

final_pattern = initial_pattern[ticks]
envelope = final_pattern.layers()[1]
initial_stator_on = final_pattern.layers()[2]
intial_stator_population = initial_stator_on.population
rotor = lt.pattern('','bs8')
rotor += envelope - initial_stator_on

# Rotor cells that have a stator cell as a neighbour
adjacent_rotor = rotor - rotor[1]

# Restrict the stator to the initial bounding  
# box with an additional 1-cell-thick border.
stator = lt.pattern('','b3s23')
stator[-1:width+1, -1:height+1] = 1
stator -= rotor

rotor_phase_0 = lt.pattern('','b3s23')
rotor_phase_0 += initial_pattern.layers()[2]
rotor_phases = [rotor_phase_0[t] & rotor for t in range(ticks + 1)]

# For each generation calculate the envelope covering all cells 
# that changed in the last tick or that border a changed cell.
# change_envelopes[0] assumes all adjacent rotor cells changed.
change_envelopes = []
for t in range(ticks):
    changed_cells = lt.pattern('','b12345678s012345678')
    changed_cells += rotor_phases[t] ^ rotor_phases[t+1]
    change_envelopes.append(changed_cells[1][-2:width+2, -2:height+2])
changed_cells = lt.pattern('','b12345678s012345678')
changed_cells += adjacent_rotor
change_envelopes.insert(0, changed_cells[1][-2:width+2, -2:height+2])

# Coordinate lists of various groups of cells
stator_cells = stator[-1:width+1, -1:height+1].coords().tolist()
stator_cells = set(map(tuple, stator_cells))
#adjacent_rotor_cells = adjacent_rotor[-1:width+1, -1:height+1].coords().tolist()
#adjacent_stator_cells = (stator & change_envelopes[0]).coords().tolist()
non_adjacent_stator_cells = (stator - change_envelopes[0]).coords().tolist()

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

# Life is the only supported rule at the moment.
rule_is_life = True

print("\nSetting up optimization search...")

model = cp_model.CpModel()

stator_int = {}
stator_bool = {}
for x,y in stator_cells:
    stator_int[x,y] = model.NewIntVar(0,1,f'stator_int_{x}_{y}')
    stator_bool[x,y] = model.NewBoolVar(f'stator_bool_{x}_{y}')
    model.Add(stator_int[x,y] == 1).OnlyEnforceIf(stator_bool[x,y])
    model.Add(stator_int[x,y] == 0).OnlyEnforceIf(stator_bool[x,y].Not())

    # Boundary stator cells must be OFF
    if x == -1 or x == width or y == -1 or y == height:
        model.Add(stator_int[x,y] == 0)

# Apply the CA rules to the stator cells
for x,y,rotor_sum in stator_neighbour_counts:
    stator_sum = sum(   stator_int[x+u, y+v]
                        for u in [-1,0,1] for v in [-1,0,1]
                        if (u,v) != (0,0) and (x+u, y+v) in stator_cells
                    )
    neighbour_sum = rotor_sum + stator_sum
    if rule_is_life:
        model.AddLinearConstraint(neighbour_sum, 2, 3).OnlyEnforceIf(stator_bool[x,y])
        model.Add(neighbour_sum != 3).OnlyEnforceIf(stator_bool[x,y].Not())
    else:
        for count in range(9):
            if B[count]:
                model.Add(neighbour_sum != count).OnlyEnforceIf(stator_bool[x,y].Not())
            if not S[count]:
                model.Add(neighbour_sum != count).OnlyEnforceIf(stator_bool[x,y])

# Apply the CA rules to the rotor transitions
for x,y,rotor_sum,state_0,state_1 in rotor_transitions:
    stator_sum = sum(   stator_int[x+u, y+v]
                        for u in [-1,0,1] for v in [-1,0,1]
                        if (u,v) != (0,0) and (x+u, y+v) in stator_cells
                    )
    neighbour_sum = rotor_sum + stator_sum
    if rule_is_life:
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
            if state_0 == 0 and state_1 == 0 and B[count]:
                model.Add(neighbour_sum != count)
            if state_0 == 0 and state_1 == 1 and not B[count]:
                model.Add(neighbour_sum != count)
            if state_0 == 1 and state_1 == 0 and S[count]:
                model.Add(neighbour_sum != count)
            if state_0 == 1 and state_1 == 1 and not S[count]:
                model.Add(neighbour_sum != count)

# Use the input pattern as a hint to the solver
intial_stator_on_cells = set(map(tuple, (stator & initial_stator_on).coords().tolist()))
intial_stator_off_cells = set(map(tuple, (stator - initial_stator_on).coords().tolist()))
for x,y in intial_stator_on_cells:
    model.AddHint(stator_int[x,y], 1)
for x,y in intial_stator_off_cells:
    model.AddHint(stator_int[x,y], 1)

size = sum(stator_int[x,y] for x,y in stator_cells)
#model.Add(size <= intial_stator_population)
model.Minimize(size)

solver = cp_model.CpSolver()

print("Beginning optimization search...")

status = solver.Solve(model)

if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
    if status == cp_model.OPTIMAL:
        print("\nOptimal solution found.")
    elif status == cp_model.FEASIBLE:
        print("\nNon-optimal solution found.")

    output_pattern = lt.pattern('','b3s23')
    output_pattern += initial_pattern & rotor

    for x,y in stator_cells:
        if solver.Value(stator_int[x,y]) == 1:
            output_pattern[x,y] = 1
    if output_pattern.population < initial_population:
        print(f"\nFinal population: {output_pattern.population}")
        print(output_pattern.rle_string())
    elif status == cp_model.OPTIMAL:
        print("\nInput pattern is already optimal.")
    else:
        print("\nNew pattern is not smaller than the input pattern.")
elif status == cp_model.UNKNOWN:
    print("\nUNKNOWN: no solutions were found, but the model was not proven infeasible.")
elif status == cp_model.INFEASIBLE:
    print(status)
    print("\nINFEASIBLE: no solutions possible.")
