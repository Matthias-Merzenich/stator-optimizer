#!/usr/bin/env python3

from ortools.sat.python import cp_model

from stator_opt.analysis import get_transition_sets
from stator_opt.environment import validate_lifelib
from stator_opt.optimizer import apply_conditional_transform, build_model
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

    verbose_print(f"Initial population: {the_pattern.initial_population}")
    verbose_print(f"Initial bounding box:"
                  f" {the_pattern.width} x {the_pattern.height}")
    verbose_print(f"Search bounding box:"
                  f" {the_pattern.width + adjust_left + adjust_right}"
                  f" x {the_pattern.height + adjust_top + adjust_bottom}\n")

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
    boundary_cells = the_pattern.stator_boundary.coords().tolist()
    boundary_cells = set(map(tuple, boundary_cells))

    stator_neighbor_counts, rotor_transitions = get_transition_sets(the_pattern)

    verbose_print("Setting up optimization search...")
    model = cp_model.CpModel()
    stator_int = {}
    stator_bool = {}
    build_model(
        rulestring,
        model,
        stator_int,
        stator_bool,
        stator_cells,
        boundary_cells,
        args.boundary,
        forced_on_cells,
        forced_off_cells,
        stator_neighbor_counts,
        rotor_transitions
    )

    # Transformation sets that generate the given symmetry
    FORCED_TRANSFORMATIONS = {
        "C1": [],
        "C2": ["rotate_180"],
        "C4": ["rotate_90"],
        "D2-": ["flip_y"],
        "D2|": ["flip_x"],
        "D2/": ["flip_diag"],
        "D2\\": ["flip_reverse_diag"],
        "D4+": ["flip_x", "flip_y"],
        "D4X": ["flip_diag", "flip_reverse_diag"],
        "D8": ["flip_x", "flip_diag"]
    }

    trans_bool = {}
    symm_bool = {}

    TRANSFORMATION_LIST = ["rotate_90", "rotate_180", "flip_x", "flip_y",
                           "flip_diag", "flip_reverse_diag"]
    for trans in TRANSFORMATION_LIST:
        trans_bool[trans] = (model.NewBoolVar(f"trans_bool_{trans}"))

    SYMMETRY_LIST = ["C2", "C4", "D2-", "D2|", "D2/",
                     "D2\\", "D4+", "D4X", "D8"]
    for symm in SYMMETRY_LIST:
        symm_bool[symm] = model.NewBoolVar(f"symm_bool_{symm}")
        for trans in FORCED_TRANSFORMATIONS[symm]:
            model.AddImplication(symm_bool[symm], trans_bool[trans])

    transformations = list(FORCED_TRANSFORMATIONS[args.symmetry])
    if args.prefer_higher_symmetry:
        transformations = TRANSFORMATION_LIST

    for transformation in transformations:
        apply_conditional_transform(
            model,
            transformation,
            trans_bool,
            stator_int,
            stator_cells,
            stator_array
        )

    # Cause transformations to be applied unconditionally
    for trans in FORCED_TRANSFORMATIONS[args.symmetry]:
        model.Add(trans_bool[trans] == 1)

    # Use the input pattern as a hint to the solver
    intial_stator_on_cells = set(map(tuple, (   the_pattern.stator
                                                & the_pattern.initial_stator_on
                                            ).coords().tolist()))
    intial_stator_off_cells = set(map(tuple, (  the_pattern.stator
                                                - the_pattern.initial_stator_on
                                             ).coords().tolist()))
    for x,y in intial_stator_on_cells:
        model.AddHint(stator_int[x,y], 1)
    for x,y in intial_stator_off_cells:
        model.AddHint(stator_int[x,y], 1)

    stator_population = sum(stator_int[x,y] for x,y in stator_cells)
    symmetry_objective = sum(symm_bool[symm] for symm in SYMMETRY_LIST)

    model.Minimize(20 * stator_population - symmetry_objective)

    solver = cp_model.CpSolver()

    verbose_print("Beginning optimization search...\n")

    status = solver.Solve(model)

    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        if status == cp_model.OPTIMAL:
            verbose_print("Optimal solution found.\n")
        elif status == cp_model.FEASIBLE:
            verbose_print("Non-optimal solution found.\n")

        output_pattern = the_pattern.initial_two_state & the_pattern.rotor
        for x,y in stator_cells:
            if solver.Value(stator_int[x,y]) == 1:
                output_pattern[x,y] = 1

        # Print the pattern if at least one of the following conditions holds:
        # * The solution has a smaller population than the original.
        # * The original stator contains a cell outside the search area.
        # * A cell in the original stator is forced off in the solution.
        # * A cell not in the original stator is forced on in the solution.
        if (
            output_pattern.population < the_pattern.initial_population
            or (the_pattern.initial_stator_on
                - (the_pattern.stator - the_pattern.stator_boundary)).nonempty()
            or (forced_off & the_pattern.initial_stator_on).nonempty()
            or (forced_on - the_pattern.initial_stator_on).nonempty()
        ):
            verbose_print(f"Final population: {output_pattern.population}")
            print(clean_rle(output_pattern.rle_string()))
        elif status == cp_model.OPTIMAL:
            verbose_print("Input pattern is already optimal.")
        else:
            verbose_print("New pattern is not smaller than the input pattern.")
    elif status == cp_model.UNKNOWN:
        verbose_print("UNKNOWN: no solutions were found, but the model"
                      "was not proven infeasible.")
    elif status == cp_model.INFEASIBLE:
        verbose_print("INFEASIBLE: no solution possible.")
    verbose_print()


if __name__ == "__main__":
    main()
else:
    lifelib = validate_lifelib(need_latest_lifelib=True)
