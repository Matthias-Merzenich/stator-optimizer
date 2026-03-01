from .rules import get_rule


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
