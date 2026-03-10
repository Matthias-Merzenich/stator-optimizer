from .rules import get_rule
from .symmetry import symmetrize, transform


# Build the base model that applies the CA rules to
# the stator and adjacent rotor.
def build_model(rulestring,
                model,
                stator_cells,
                unforced_boundary_cells,
                forced_on_cells,
                forced_off_cells,
                stator_neighbor_counts,
                rotor_transitions,
                extra_boundary_cells,
                validate_boundary
):
    birth_at, survival_at = get_rule(rulestring)

    stator_vars = {}
    valid_boundary_vars = {}
    for x,y in stator_cells:
        stator_vars[x,y] = model.NewBoolVar(f"stator_{x}_{y}")

        # If a forced-on cell is in the forced-off boundary, we
        # prioritize the forced-off setting, because the cell
        # is outside of the search area.
        if (x,y) in forced_off_cells:
            model.Add(stator_vars[x,y] == 0)
        elif (x,y) in forced_on_cells:
            model.Add(stator_vars[x,y] == 1)

    def get_stator_sum(x, y):
        return sum( stator_vars[x+u, y+v]
                    for u in [-1,0,1] for v in [-1,0,1]
                    if (u,v) != (0,0) and (x+u, y+v) in stator_cells
                  )

    # Apply the CA rules to the stator cells
    for x,y,rotor_sum in stator_neighbor_counts:
        cell_on = [stator_vars[x,y]]
        cell_off = [stator_vars[x,y].Not()]
        if (x,y) in unforced_boundary_cells:
            if validate_boundary:
                # If valid_boundary_vars[x,y] is true,
                # then the CA rules apply at (x,y).
                valid_boundary_vars[x,y] = (
                    model.NewBoolVar(f"valid_boundary_{x}_{y}")
                )
                cell_on.append(valid_boundary_vars[x,y])
                cell_off.append(valid_boundary_vars[x,y])
            else:
                continue

        neighbor_sum = get_stator_sum(x,y) + rotor_sum
        if rulestring == "b3s23":
            model.Add(neighbor_sum != 3).OnlyEnforceIf(cell_off)
            model.AddLinearConstraint(neighbor_sum, 2, 3).OnlyEnforceIf(cell_on)
        else:
            for count in range(9):
                if birth_at[count]:
                    model.Add(neighbor_sum != count).OnlyEnforceIf(cell_off)
                if not survival_at[count]:
                    model.Add(neighbor_sum != count).OnlyEnforceIf(cell_on)

    # Apply the CA rules to the rotor transitions
    for x,y,rotor_sum,state_0,state_1 in rotor_transitions:
        neighbor_sum = get_stator_sum(x, y) + rotor_sum
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

    # If we are validating the boundary, make sure to identify when the
    # unchecked boundary would cause births outside of the search area.
    if validate_boundary:
        for x,y in extra_boundary_cells:
            valid_boundary_vars[x,y] = (
                model.NewBoolVar(f"valid_boundary_{x}_{y}")
            )
            neighbor_sum = get_stator_sum(x,y)
            for count in range(9):
                if birth_at[count]:
                    model.Add(
                        neighbor_sum != count
                    ).OnlyEnforceIf(valid_boundary_vars[x,y])

    return stator_vars, valid_boundary_vars


# For the given transformation, add conditions so that the
# pattern is invariant under the transformation if the
# variable trans_bool[transformation] is True. The actual
# state of trans_bool[transformation] is determined during
# solving or set at run time.
def apply_conditional_transform(
    model,
    transformation,
    trans_bool,
    stator_vars,
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
            stator_vars[tuple(cell_1)] == stator_vars[tuple(cell_2)]
        ).OnlyEnforceIf(trans_bool[transformation])

    stator_cells = set(stator_vars)
    symm_stator_cells = set(map(tuple, symm_stator_array.tolist()))
    for cell in stator_cells - symm_stator_cells:
        model.Add(
            stator_vars[cell] == 0
        ).OnlyEnforceIf(trans_bool[transformation])
