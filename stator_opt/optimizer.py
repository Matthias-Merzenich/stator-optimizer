from .rules import get_rule
from .symmetry import symmetrize, transform


def build_model(rulestring,
                model,
                stator_int,
                stator_bool,
                stator_cells,
                boundary_cells,
                boundary_state,
                forced_on_cells,
                forced_off_cells,
                stator_neighbor_counts,
                rotor_transitions
               ):
    birth_at, survival_at = get_rule(rulestring)

    for x,y in stator_cells:
        stator_int[x,y] = model.NewIntVar(0,1,f"stator_int_{x}_{y}")
        stator_bool[x,y] = model.NewBoolVar(f"stator_bool_{x}_{y}")
        model.Add(stator_int[x,y] == 1).OnlyEnforceIf(stator_bool[x,y])
        model.Add(stator_int[x,y] == 0).OnlyEnforceIf(stator_bool[x,y].Not())

        # If a forced-on cell is in the forced-off boundary, we
        # prioritize the forced-off setting, because the cell
        # is outside of the search area.
        if (x,y) in forced_off_cells:
            model.Add(stator_int[x,y] == 0)
        elif (x,y) in forced_on_cells:
            model.Add(stator_int[x,y] == 1)

    # Apply the CA rules to the stator cells
    for x,y,rotor_sum in stator_neighbor_counts:
        if (x,y) in boundary_cells and boundary_state == "any":
            continue
        stator_sum = sum(   stator_int[x+u, y+v]
                            for u in [-1,0,1] for v in [-1,0,1]
                            if (u,v) != (0,0) and (x+u, y+v) in stator_cells
                        )
        neighbor_sum = rotor_sum + stator_sum
        if rulestring == "b3s23":
            model.AddLinearConstraint(neighbor_sum, 2, 3).OnlyEnforceIf(
                stator_bool[x,y]
            )
            model.Add(neighbor_sum != 3).OnlyEnforceIf(
                stator_bool[x,y].Not()
            )
        else:
            for count in range(9):
                if birth_at[count]:
                    model.Add(neighbor_sum != count).OnlyEnforceIf(
                        stator_bool[x,y].Not()
                    )
                if not survival_at[count]:
                    model.Add(neighbor_sum != count).OnlyEnforceIf(
                        stator_bool[x,y]
                    )

    # Apply the CA rules to the rotor transitions
    for x,y,rotor_sum,state_0,state_1 in rotor_transitions:
        stator_sum = sum(   stator_int[x+u, y+v]
                            for u in [-1,0,1] for v in [-1,0,1]
                            if (u,v) != (0,0) and (x+u, y+v) in stator_cells
                        )
        neighbor_sum = rotor_sum + stator_sum
        if rulestring == "b3s23":
            if state_0 == 0 and state_1 == 0:
                model.Add(neighbor_sum != 3)
            if state_0 == 0 and state_1 == 1:
                model.Add(neighbor_sum == 3)
            if state_0 == 1 and state_1 == 0:
                model.Add(neighbor_sum != 2)
                model.Add(neighbor_sum != 3)
            if state_0 == 1 and state_1 == 1:
                model.AddLinearConstraint(neighbor_sum, 2, 3)
        else:
            for count in range(9):
                if state_0 == 0 and state_1 == 0 and birth_at[count]:
                    model.Add(neighbor_sum != count)
                if state_0 == 0 and state_1 == 1 and not birth_at[count]:
                    model.Add(neighbor_sum != count)
                if state_0 == 1 and state_1 == 0 and survival_at[count]:
                    model.Add(neighbor_sum != count)
                if state_0 == 1 and state_1 == 1 and not survival_at[count]:
                    model.Add(neighbor_sum != count)


def apply_conditional_transform(
    model,
    transformation,
    trans_bool,
    stator_int,
    stator_cells,
    stator_array
):
    SYMM_FROM_TRANS = {
        "flip_x": "D2|",
        "flip_y": "D2-",
        "flip_diag": "D2/",
        "flip_reverse_diag": "D2\\",
        "rotate_90": "C4",
        "rotate_180": "C2"
    }
    symm_stator_array, center = symmetrize(
        stator_array, SYMM_FROM_TRANS[transformation]
    )
    transformed_coords = transform( symm_stator_array,
                                    transformation,
                                    center
                                  ).tolist()
    for cell_1, cell_2 in zip(symm_stator_array.tolist(), transformed_coords):
        model.Add(
            stator_int[tuple(cell_1)] == stator_int[tuple(cell_2)]
        ).OnlyEnforceIf(trans_bool[transformation])

    symm_stator_cells = set(map(tuple, symm_stator_array.tolist()))
    for cell in stator_cells - symm_stator_cells:
        model.Add(
            stator_int[cell] == 0
        ).OnlyEnforceIf(trans_bool[transformation])
