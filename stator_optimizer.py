#!/usr/bin/env python3

from ortools.sat.python import cp_model

from stator_opt.analysis import get_transition_sets
from stator_opt.environment import validate_lifelib
from stator_opt.objectives import (
    SYMMETRY_WEIGHTS,
    ObjectiveStats,
    apply_box_objectives,
    apply_symmetry_objective,
    get_objective_expression,
)
from stator_opt.optimizer import build_model
from stator_opt.parsers import (
    parse_arguments,
    parse_rotor_descriptor,
    read_file,
)
from stator_opt.pattern_stats import PatternStats
from stator_opt.rules import build_ruletable, get_rule, get_rulestring
from stator_opt.setup import create_ruletable_session, get_pattern

# Before importing lifelib, we must run some checks to see if the
# features used are compatible with the current lifelib version.
# If not, validate_lifelib() downloads the latest version from
# https://gitlab.com/apgoucher/lifelib
# Specifically, if the input file is in rotor descriptor format, we
# need lifelib version 2.5.9 or higher.


def clean_rle(rle):
    return "\n".join([line for line in rle.splitlines()
                      if not line.startswith('#')])


def main():
    args = parse_arguments()
    adjust_left, adjust_right, adjust_top, adjust_bottom = args.adjust
    verbose_print = print if not args.solution_only else lambda *a, **k: None

    unique_objectives = []
    for objective in args.objectives:
        if objective not in unique_objectives:
            unique_objectives.append(objective)
    if "max_pop" in unique_objectives and "min_pop" in unique_objectives:
        max_pop_index = unique_objectives.index("max_pop")
        min_pop_index = unique_objectives.index("min_pop")
        unique_objectives.pop(max(max_pop_index, min_pop_index))
    args.objectives = unique_objectives

    rle = read_file(args.input_file)

    rulestring, grid = parse_rotor_descriptor(rle)
    if grid is not None:
        lifelib = validate_lifelib(need_latest_lifelib=True)

        rulestring = lifelib.sanirule(rulestring)
        get_rule(rulestring)    # Validate the rule.
        ruletable = build_ruletable(rulestring)

        sess = create_ruletable_session(rulestring, ruletable, lifelib)
        lt = sess.lifetree(n_layers=1)
        lt5 = sess.lifetree(n_layers=5)

        grid_rle = '$'.join(grid)
        grid_rle = grid_rle.translate(str.maketrans("012345678@ABCDEFGH",
                                                    "BDFHJLNPRACEGIKMOQ"
                                                   ))
        grid_rle = f"x = 0, y = 0, rule=Rotor{rulestring}\n{grid_rle}!"
        initial_pattern = lt5.pattern(grid_rle)
        forced_on = lt.pattern("", rulestring)
        forced_off = lt.pattern("", rulestring)
    else:
        lifelib = validate_lifelib()

        rulestring = get_rulestring(rle, lifelib.load_rules("life").lifetree())
        get_rule(rulestring)    # Validate the rule.

        sess = lifelib.load_rules(
            rulestring,
            "bs8",
            "b12345678s012345678",
            "b1e2-cn3-c4-c5678s012345678"
        )
        lt = sess.lifetree(n_layers=1)
        lt5 = sess.lifetree(n_layers=5)

        initial_pattern, forced_on, forced_off = get_pattern(rle, lt, lt5)

    the_pattern = PatternStats(initial_pattern, lt)

    verbose_print("Analyzing pattern...")

    the_pattern.analyze_pattern(
        args.ticks,
        args.adjust,
        args.distance
    )

    forced_on &= the_pattern.stator
    forced_off &= the_pattern.stator

    # For still life searches, the distance argument is
    # applied to the forced cells, rather than the rotor.
    if args.distance is not None and the_pattern.rotor.empty():
        final_mask = the_pattern.stator - the_pattern.stator_boundary
        the_pattern.stator = lt.pattern("", "b1e2-cn3-c4-c5678s012345678")
        the_pattern.stator += forced_on + forced_off
        the_pattern.stator = the_pattern.stator[args.distance]
        the_pattern.stator &= final_mask

        the_pattern.stator_boundary = lt.pattern("", "b12345678s012345678")
        the_pattern.stator_boundary += the_pattern.stator
        the_pattern.stator_boundary = (
            the_pattern.stator_boundary[1] - the_pattern.stator
        )
        the_pattern.stator += the_pattern.stator_boundary

    if args.boundary == "off":
        forced_off += the_pattern.stator_boundary

    forced_on_cells = forced_on.coords().tolist()
    forced_on_cells = set(map(tuple, forced_on_cells))
    forced_off_cells = forced_off.coords().tolist()
    forced_off_cells = set(map(tuple, forced_off_cells))
    stator_array = the_pattern.stator.coords()
    stator_cells = stator_array.tolist()
    stator_cells = set(map(tuple, stator_cells))
    if not stator_cells:
        raise RuntimeError("Search area is empty.")
    boundary_cells = the_pattern.stator_boundary.coords().tolist()
    boundary_cells = set(map(tuple, boundary_cells))

    objective_stats = ObjectiveStats(the_pattern)
    verbose_print(objective_stats.get_stats_string(args.objectives))

    stator_neighbor_counts, rotor_transitions = get_transition_sets(the_pattern)

    verbose_print("Setting up optimization search...")
    model = cp_model.CpModel()

    # Apply the life rules for the stator cells and rotor transitions.
    stator_vars = build_model(
        rulestring,
        model,
        stator_cells,
        boundary_cells,
        args.boundary,
        forced_on_cells,
        forced_off_cells,
        stator_neighbor_counts,
        rotor_transitions
    )

    # Apply constraints relating to the bounding box and bounding diamond.
    box_vars, objective_dict = apply_box_objectives(
        model, the_pattern, args.objectives, stator_vars, stator_array
    )

    # Apply constraints relating to symmetry.
    symm_bool = apply_symmetry_objective(
        model, args.objectives, args.symmetry, stator_vars, stator_array
    )

    # Finish building objective_dict, which is used to construct the
    # expression that we want to minimize.
    max_pop = len(stator_cells)
    max_symm = sum(SYMMETRY_WEIGHTS.values())
    pop_var = model.NewIntVar(0, max_pop, 'stator_pop')
    symm_var = model.NewIntVar(0, max_symm, 'symmetry_sum')
    model.Add(pop_var == cp_model.LinearExpr.Sum(list(stator_vars.values())))
    model.Add(symm_var == cp_model.LinearExpr.WeightedSum(
        [symm_bool[symm] for symm in SYMMETRY_WEIGHTS],
        [SYMMETRY_WEIGHTS[symm] for symm in SYMMETRY_WEIGHTS]
    ))
    objective_dict.update({
        "min_pop":  [pop_var, "min", max_pop],
        "max_pop":  [pop_var, "max", max_pop],
        "symmetry": [symm_var, "max", max_symm]
    })

    total_objective, init_objective = get_objective_expression(
        the_pattern, objective_stats, objective_dict, args.objectives
    )
    model.Minimize(total_objective)

    solver = cp_model.CpSolver()

    verbose_print(f"  Undetermined cells:"
                  f" {len(stator_cells - forced_on_cells - forced_off_cells)}")
    verbose_print(f"  Variables: {len(model.Proto().variables)}")
    verbose_print(f"  Constraints: {len(model.Proto().constraints)}")

    verbose_print("Beginning optimization search...\n")

    status = solver.Solve(model)

    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        if status == cp_model.OPTIMAL:
            verbose_print("Optimal solution found.")
        elif status == cp_model.FEASIBLE:
            verbose_print("Non-optimal solution found.")

        # Print the pattern if at least one of the following conditions holds:
        # * The solution has a better objective value than the original.
        # * The original stator contains a cell outside the search area.
        # * A cell in the original stator is forced off in the solution.
        # * A cell not in the original stator is forced on in the solution.
        if (
            round(solver.ObjectiveValue()) < init_objective
            or (the_pattern.initial_stator_on
                - (the_pattern.stator - the_pattern.stator_boundary)).nonempty()
            or (forced_off & the_pattern.initial_stator_on).nonempty()
            or (forced_on - the_pattern.initial_stator_on).nonempty()
        ):
            output_pattern = the_pattern.initial_two_state & the_pattern.rotor
            for x,y in stator_cells:
                if solver.Value(stator_vars[x,y]) == 1:
                    output_pattern[x,y] = 1

            objective_stats.store_final_stats(
                solver, stator_vars, box_vars, symm_bool
            )
            verbose_print(objective_stats.get_stats_string(args.objectives))
            print(clean_rle(output_pattern.rle_string()))
        elif status == cp_model.OPTIMAL:
            verbose_print("\nInput pattern is already optimal.")
        else:
            verbose_print("\nNew pattern is not more optimal than the"
                          " input pattern.")
    elif status == cp_model.UNKNOWN:
        verbose_print("UNKNOWN: no solutions were found, but the model"
                      " was not proven infeasible.")
    elif status == cp_model.INFEASIBLE:
        verbose_print("INFEASIBLE: no solution possible.")
    verbose_print()


if __name__ == "__main__":
    main()
else:
    lifelib = validate_lifelib(need_latest_lifelib=True)
