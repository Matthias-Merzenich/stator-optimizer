import numpy as np
from scipy.spatial import cKDTree


def get_transition_sets(pattern):
    ticks = len(pattern.rotor_phases) - 1
    non_adjacent_stator_cells = (pattern.stator
                                 - pattern.adjacent_rotor[1]).coords().tolist()

    # Set of tuples (x,y,c) where (x,y) is a stator cell and c is a count of
    # its on rotor neighbors in some generation. A single stator cell will
    # correspond to multiple elements in this set.
    stator_neighbor_counts = set()

    # Set of tuples (x, y, c, state_0, state_1) where (x,y) is a rotor cell,
    # c is a count of its on rotor neighbors in some generation, state_0 is
    # the state of the cell in that generation, and state_1 is the state in
    # the next generation.
    rotor_transitions = set()

    for t in range(ticks):
        tree = cKDTree(pattern.rotor_phases[t].coords())

        gen_t_stator_coords = (pattern.change_envelopes[t]
                               & pattern.stator).coords()
        gen_t_stator_neighbor_counts = (
            tree.query_ball_point(
                gen_t_stator_coords, r=1, p=np.inf, return_length=True
            )
        )
        gen_t_stator_neighbor_counts = np.expand_dims(
            gen_t_stator_neighbor_counts, axis=1
        )

        coords_and_counts = np.hstack(( gen_t_stator_coords,
                                        gen_t_stator_neighbor_counts
                                      )).astype(np.int64).tolist()
        stator_neighbor_counts.update(set(map(tuple, coords_and_counts)))

        gen_t_adjacent_rotor_coords = (pattern.change_envelopes[t]
                                       & pattern.adjacent_rotor).coords()
        gen_t_rotor_state_0 = pattern.rotor_phases[t][gen_t_adjacent_rotor_coords]
        gen_t_rotor_state_1 = pattern.rotor_phases[t+1][gen_t_adjacent_rotor_coords]
        gen_t_rotor_neighbor_counts = (
            tree.query_ball_point(
                gen_t_adjacent_rotor_coords, r=1, p=np.inf, return_length=True
            )
            - tree.query_ball_point(
                gen_t_adjacent_rotor_coords, r=0, p=np.inf, return_length=True
            )
        )
        gen_t_rotor_neighbor_counts = (
            np.expand_dims(gen_t_rotor_neighbor_counts, axis=1)
        )
        gen_t_rotor_state_0 = np.expand_dims(gen_t_rotor_state_0, axis=1)
        gen_t_rotor_state_1 = np.expand_dims(gen_t_rotor_state_1, axis=1)

        coords_and_counts = np.hstack(( gen_t_adjacent_rotor_coords,
                                        gen_t_rotor_neighbor_counts,
                                        gen_t_rotor_state_0,
                                        gen_t_rotor_state_1
                                      )).astype(np.int64).tolist()
        rotor_transitions.update(set(map(tuple, coords_and_counts)))

    for x,y in non_adjacent_stator_cells:
        stator_neighbor_counts.add((x,y,0))

    return stator_neighbor_counts, rotor_transitions
