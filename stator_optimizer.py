#!/usr/bin/env python3

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from importlib.metadata import PackageNotFoundError, version

import numpy as np
from ortools.sat.python import cp_model
from scipy.spatial import cKDTree

# Before importing lifelib, we must run some checks to see if the
# features used are compatible with the current lifelib version.
# If not, the script attempts to download the latest version from
# https://gitlab.com/apgoucher/lifelib
# Specifically, if the input file is in rotor descriptor format, we
# need lifelib version 2.5.9 or higher.

def parse_arguments():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="A program to optimize the stator of patterns in Life-like cellular automata.",
        add_help=False
    )
    parser.add_argument(
        "input_file",
        help="File containing the pattern to be optimized. When the\n"
             "pattern is saved using a history rule, you can specify\n"
             "stator cell states in the solution as follows:\n"
             "------------------------------------------------\n"
             " history state | input state | solution state\n"
             "---------------+-------------+------------------\n"
             " 0 (black)     | off         | either on or off\n"
             " 1 (green)     | on          | either on or off\n"
             " 2 (blue)      | off         | off\n"
             " 3 (white)     | on          | off\n"
             " 4 (red)       | off         | on\n"
             " 5 (yellow)    | on          | on\n"
             " 6 (gray)      | off         | either on or off\n"
             "This only applies to stator cells. If a cell is\n"
             "determined to be part of the rotor, then its state\n"
             "will not be changed in the solution.\n"
             "Accepted input formats are RLE (.rle), macrocell (.mc),\n"
             "and rotor descriptor format.\n\n"
    )
    parser.add_argument(
        "ticks",
        type=int,
        help="Number of time steps to run the pattern for analysis."
    )
    parser.add_argument(
        "-h", "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit."
    )
    parser.add_argument(
        "-a", "--adjust",
        type=int,
        nargs=4,
        default=[0, 0, 0, 0],
        metavar=("LEFT", "RIGHT", "TOP", "BOTTOM"),
        help="Expand the search box by the given distances. Negative\n"
             "values contract the search box."
    )
    parser.add_argument(
        "--solution_only",
        action="store_true",
        help="Only print the solution to the optimization problem.\n"
             "If there is no solution or if the input pattern is\n"
             "already optimal, print nothing."
    )
    return parser.parse_args()


def read_input_file(file_name):
    try:
        with open(file_name, 'r') as file:
            rle = file.read()
    except FileNotFoundError:
        raise FileNotFoundError("Pattern file not found.")

    return rle


def parse_version(version_string):
    return tuple(map(int, (version_string.split('.') + ['0', '0'])[:3]))


def validate_lifelib():
    min_version = "2.5.9"
    repo_url = "https://gitlab.com/apgoucher/lifelib.git"
    clone_required = False
    try:
        current_version = version("python-lifelib")
        if (
            parse_version(current_version) < parse_version(min_version)
            and not os.path.exists("lifelib")
        ):
            print("Incompatible lifelib version: "
                  f"{current_version} < {min_version}", file=sys.stderr
                 )
            clone_required = True
    except PackageNotFoundError:
        if not os.path.exists("lifelib"):
            print("lifelib not found.", file=sys.stderr)
            clone_required = True

    if clone_required:
        print("Fetching latest version...", file=sys.stderr)
        if shutil.which("git") is None:
                print("ERROR: 'git' is not installed or is not in PATH.", file=sys.stderr)
                print("Please install git or download the "
                      "latest version of lifelib from\n"
                      "https://gitlab.com/apgoucher/lifelib", file=sys.stderr
                     )
        try:
            subprocess.run(["git", "clone", repo_url, "lifelib"], check=True)
            print("Clone successful.", file=sys.stderr)
        except subprocess.CalledProcessError as e:
            print(f"Error cloning repository: {e}", file=sys.stderr)
            return

    if os.path.exists("lifelib"):
        sys.path.insert(0, os.path.abspath("lifelib"))


# Returns the rulestring and rotor grid if the input string can
# be read as a rotor descriptor. Otherwise returns None, None.
def parse_rotor_descriptor(input_string):
    RE_PATTERN = re.compile(r"""
        (?:R(\S+)\s)?               # Optional Rule
        (?:p[0-9]+\s)?              # Optional period (ignored)
        (?:r[0-9]+\s)?              # Optional rotor size (ignored)
        (?:([0-9]+)x([0-9]+)\s)?    # Optional height and width
        (                           # Rotor grid
            (?:
                (?![ \t]{2})        # Stop if 2+ horizontal whitespace chars.
                [ \n\r.0-8@A-H]     # Allowed grid characters: ., 0-8, @, A-H
            )+
        )
    """, re.VERBOSE)
    match = RE_PATTERN.match(input_string)
    if not match:
        return None, None
    rulestring, height, width, grid = match.groups()
    grid = grid.split()
    if not grid:
        return None, None
    if height is None:
        height = len(grid)
        width = len(grid[0])
    else:
        height = int(height)
        width = int(width)
    if len(grid) < height or any(len(row) != width for row in grid[:height]):
        raise ValueError("Rotor descriptor dimensions do not match grid.")
    if set(''.join(grid[:height])) == {"."}:
        return None, None
    if rulestring is None:
        rulestring = "B3/S23"

    return rulestring, grid


if __name__ == "__main__":
    if parse_rotor_descriptor(read_input_file(parse_arguments().input_file))[0] is not None:
        validate_lifelib()
else:
    validate_lifelib()

import lifelib


class PatternStats:
    def __init__(self, initial_pattern, lt):
        self.initial = initial_pattern
        _, _, self.width, self.height = self.initial.bounding_box
        rulestring = get_rulestring_from_pattern(self.initial)
        if rulestring.startswith("xrotor"):
            rulestring = rulestring[len("xrotor"):]
        else:
            # Remove state-2 cells
            self.initial = ((self.initial >> 2) << 2) + (self.initial >> 2)
        self.initial_two_state = lt.pattern("", rulestring)
        self.initial_two_state += self.initial
        self.rotor = lt.pattern("", "bs8")
        self.stator = lt.pattern("", "b12345678s012345678")
        self.initial_stator_on = lt.pattern("", "b12345678s012345678")
        self.stator_boundary = None     # The boundary of the search area.
        self.adjacent_rotor = None      # Rotor cells with a stator neighbour.
        self.rotor_phases = None
        self.change_envelopes = None

    @property
    def initial_population(self):
        if self.initial.getrule().endswith("History"):
            return self.initial_two_state.population
        else:
            return np.inf

    def make_stator(self, adjustments):
        adjust_left, adjust_right, adjust_top, adjust_bottom = adjustments

        # Restrict the stator to the search box
        # with an additional 1-cell-thick border.
        self.stator[    0 - adjust_left : self.width + adjust_right,
                        0 - adjust_top : self.height + adjust_bottom
                   ] = 1
        self.stator -= self.rotor
        self.stator_boundary = self.stator[1] - self.stator - self.rotor
        self.stator += self.stator_boundary

        # Note that self.stator[1] comes first in the assignment, because
        # we want self.adjacent_rotor to have rule B12345678/S012345678.
        self.adjacent_rotor = self.stator[1] & (self.rotor - self.rotor[1])

    def analyze_pattern(self, ticks, adjustments):
        pattern_is_history = self.initial.getrule().endswith("History")
        if pattern_is_history:
            final_pattern = self.initial[ticks]
            envelope = flatten_multistate(final_pattern)
            self.rotor += flatten_multistate(final_pattern ^ self.initial)
            self.initial_stator_on += envelope ^ self.rotor
        else:
            self.rotor += flatten_multistate(self.initial)

        self.make_stator(adjustments)

        previous_phase = self.initial
        rotor_mask = self.adjacent_rotor[1] & self.rotor
        self.rotor_phases = [rotor_mask & self.initial]
        for t in range(1, ticks+1):
            this_phase = previous_phase[1]
            self.rotor_phases.append(rotor_mask & this_phase)
            previous_phase = this_phase

        # For each generation calculate the envelope covering all cells
        # that changed in the last tick or that border a changed cell.
        # change_envelopes[0] assumes all adjacent rotor cells changed.
        rotor_mask = self.adjacent_rotor[1]
        self.change_envelopes = [rotor_mask]
        for t in range(ticks):
            changed_cells = self.rotor_phases[t] ^ self.rotor_phases[t+1]
            self.change_envelopes.append(changed_cells[1] & rotor_mask)


def flatten_multistate(pattern):
    n = pattern.owner.n_layers
    flattened_pattern = pattern
    for shift in range(n):
        flattened_pattern = flattened_pattern + (pattern >> shift)
    final_shift = 2**((n-1).bit_length()) - 1
    return (flattened_pattern << final_shift) >> final_shift


def clean_rle(rle):
    return "\n".join([line for line in rle.splitlines() if not line.startswith('#')])


def get_rulestring_from_pattern(pattern):
    rulestring = pattern.getrule()
    if rulestring.endswith("History"):
        rulestring = rulestring[:-len("History")]
    return rulestring


def get_rulestring(rle):
    return get_rulestring_from_pattern(lifelib.load_rules("b3s23").lifetree().pattern(rle))


def get_rule(rulestring):
    match = re.match(r"b([2-8]*)s([0-8]*)$", rulestring)
    if not match:
        raise ValueError(f"Only Life-like rules without B0 or B1 are supported. Got: {rulestring}")
    birth_digits, survival_digits = match.groups()

    if not survival_digits:
        raise ValueError("Invalid rule: no survival conditions.")

    def increasing(s): return all(s[i] < s[i+1] for i in range(len(s)-1))
    if not increasing(birth_digits) or not increasing(survival_digits):
        raise ValueError(f"Only Life-like rules are supported. Got: {rulestring}")

    birth_at = [False] * 9
    survival_at = [False] * 9
    for digit in birth_digits:
        birth_at[int(digit)] = True
    for digit in survival_digits:
        survival_at[int(digit)] = True

    return birth_at, survival_at


def build_ruletable(rulestring):
    birth_at, survival_at = get_rule(rulestring)
    ruletable = (f"@RULE Rotor{rulestring}\n"
                 "@TABLE\n"
                 "n_states:19\n"
                 "neighborhood:Moore\n"
                 "symmetries:permute\n"
                )

    for i in range(8):
        ruletable += f"var on{i} = {{1,3,5,7,9,11,13,15,17}}\n"
        ruletable += f"var off{i} = {{0,2,4,6,8,10,12,14,16,18}}\n"

    def stator_neighbours(state):
        return (state - 1) // 2

    for total_neighbours in range(9):
        for state in range(1,19):
            rotor_neighbours = total_neighbours - stator_neighbours(state)
            if rotor_neighbours < 0:
                continue
            if state % 2 == 1 and not survival_at[total_neighbours]:
                ruletable += f"{state}, "
                for i in range(8):
                    if i < rotor_neighbours:
                        ruletable += f"on{i},"
                    else:
                        ruletable += f"off{i},"
                ruletable += f" {state + 1}\n"
            if state % 2 == 0 and birth_at[total_neighbours]:
                ruletable += f"{state} "
                for i in range(8):
                    if i < rotor_neighbours:
                        ruletable += f"on{i},"
                    else:
                        ruletable += f"off{i},"
                ruletable += f" {state - 1}\n"

    return ruletable


def create_ruletable_session(rulestring, ruletable):
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as ruletable_file:
            ruletable_file.write(ruletable)
            ruletable_path = ruletable_file.name
    try:
        sess = lifelib.load_rules(rulestring, "bs8", "b12345678s012345678", ruletable_path)
    finally:
        if os.path.exists(ruletable_path):
            os.remove(ruletable_path)

    return sess


def get_pattern(rle, lt, lt5):
    rulestring = get_rulestring(rle)

    # Convert the input file to a history RLE. It's faster to change the
    # cell states by modifying the RLE than to change them in the lifetree.
    initial_pattern = lt5.pattern("", rulestring + "History")
    initial_pattern += lt5.pattern(rle)
    rle = initial_pattern.rle_string()

    # Identify cells to be forced on or off in the solution
    rle_forced_on = rle.translate(str.maketrans("ABCF" + "DE",
                                                "BBBB" + "AA"
                                               ))
    rle_forced_off = rle.translate(str.maketrans("ADEF" + "BC",
                                                 "BBBB" + "AA"
                                                ))
    forced_on = lt.pattern("", rulestring) + lt5.pattern(rle_forced_on)
    forced_off = lt.pattern("", rulestring) + lt5.pattern(rle_forced_off)

    # Convert all on cells to state 5 and all off cells to state 2.
    # In LifeHistory, State 5 cells will only remain state 5 as long as
    # they never die, so they can be used to detect the initial stator.
    rle = rle.translate(str.maketrans("BDF" + "ACE",
                                      "BBB" + "EEE"
                                     ))
    initial_pattern = lt5.pattern(rle)
    if initial_pattern.empty():
        raise ValueError("Input pattern is empty.")

    return initial_pattern, forced_on, forced_off


def get_transition_sets(pattern):
    ticks = len(pattern.rotor_phases) - 1
    non_adjacent_stator_cells = (pattern.stator - pattern.adjacent_rotor[1]).coords().tolist()

    # Set of tuples (x,y,c) where (x,y) is a stator cell and c is a count of
    # its on rotor neighbours in some generation. A single stator cell will
    # correspond to multiple elements in this set.
    stator_neighbour_counts = set()

    # Set of tuples (x, y, c, state_0, state_1) where (x,y) is a rotor cell,
    # c is a count of its on rotor neighbours in some generation, state_0 is
    # the state of the cell in that generation, and state_1 is the state in
    # the next generation.
    rotor_transitions = set()

    for t in range(ticks):
        tree = cKDTree(pattern.rotor_phases[t].coords())

        gen_t_stator_coords = (pattern.change_envelopes[t] & pattern.stator).coords()
        gen_t_stator_neighbour_counts = (
            tree.query_ball_point(gen_t_stator_coords, 1, p=np.inf, return_length=True)
        )
        gen_t_stator_neighbour_counts = np.expand_dims(gen_t_stator_neighbour_counts, axis=1)

        coords_and_counts = np.hstack(( gen_t_stator_coords,
                                        gen_t_stator_neighbour_counts
                                      )).astype(np.int64).tolist()
        stator_neighbour_counts.update(set(map(tuple, coords_and_counts)))

        gen_t_adjacent_rotor_coords = (pattern.change_envelopes[t] & pattern.adjacent_rotor).coords()
        gen_t_rotor_state_0 = pattern.rotor_phases[t][gen_t_adjacent_rotor_coords]
        gen_t_rotor_state_1 = pattern.rotor_phases[t+1][gen_t_adjacent_rotor_coords]
        gen_t_rotor_neighbour_counts = (
            tree.query_ball_point(gen_t_adjacent_rotor_coords, 1, p=np.inf, return_length=True)
            - tree.query_ball_point(gen_t_adjacent_rotor_coords, 0, p=np.inf, return_length=True)
        )
        gen_t_rotor_neighbour_counts = np.expand_dims(gen_t_rotor_neighbour_counts, axis=1)
        gen_t_rotor_state_0 = np.expand_dims(gen_t_rotor_state_0, axis=1)
        gen_t_rotor_state_1 = np.expand_dims(gen_t_rotor_state_1, axis=1)

        coords_and_counts = np.hstack(( gen_t_adjacent_rotor_coords,
                                        gen_t_rotor_neighbour_counts,
                                        gen_t_rotor_state_0,
                                        gen_t_rotor_state_1
                                      )).astype(np.int64).tolist()
        rotor_transitions.update(set(map(tuple, coords_and_counts)))

    for x,y in non_adjacent_stator_cells:
        stator_neighbour_counts.add((x,y,0))

    return stator_neighbour_counts, rotor_transitions


def build_model(rulestring,
                model,
                stator_int,
                stator_bool,
                stator_cells,
                forced_on_cells,
                forced_off_cells,
                stator_neighbour_counts,
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


def main(args):
    adjust_left, adjust_right, adjust_top, adjust_bottom = args.adjust
    verbose_print = print if not args.solution_only else lambda *a, **k: None

    rle = read_input_file(args.input_file)

    rulestring, grid = parse_rotor_descriptor(rle)
    if grid is not None:
        rulestring = get_rulestring(f"x = 0, y = 0, rule = {rulestring}\nb!")
        get_rule(rulestring)    # Validate the rule.
        current_version = lifelib.__version__
        if parse_version(current_version) < parse_version("2.5.9"):
            print("Failed to import latest lifelib package.", file=sys.stderr)
            sys.exit(f"Incompatible lifelib version: {current_version} < 2.5.9.\n"
                     "Rotor descriptor searches disabled.")

        ruletable = build_ruletable(rulestring)

        sess = create_ruletable_session(rulestring, ruletable)
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
        rulestring = get_rulestring(rle)
        get_rule(rulestring)    # Validate the rule.

        sess = lifelib.load_rules(rulestring, "bs8", "b12345678s012345678")
        lt = sess.lifetree(n_layers=1)
        lt5 = sess.lifetree(n_layers=5)

        initial_pattern, forced_on, forced_off = get_pattern(rle, lt, lt5)

    the_pattern = PatternStats(initial_pattern, lt)

    verbose_print(f"Initial population: {the_pattern.initial_population}")
    verbose_print(f"Initial bounding box: {the_pattern.width} x {the_pattern.height}")
    verbose_print(f"Search bounding box:"
                  f" {the_pattern.width + adjust_left + adjust_right}"
                  f" x {the_pattern.height + adjust_top + adjust_bottom}\n")

    verbose_print("Analyzing pattern...")

    the_pattern.analyze_pattern(args.ticks, args.adjust)

    forced_on &= the_pattern.stator
    forced_off &= the_pattern.stator
    forced_on_cells = forced_on.coords().tolist()
    forced_on_cells = set(map(tuple, forced_on_cells))
    forced_off_cells = (forced_off + the_pattern.stator_boundary).coords().tolist()
    forced_off_cells = set(map(tuple, forced_off_cells))
    stator_cells = the_pattern.stator.coords().tolist()
    stator_cells = set(map(tuple, stator_cells))

    stator_neighbour_counts, rotor_transitions = get_transition_sets(the_pattern)

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
        forced_on_cells,
        forced_off_cells,
        stator_neighbour_counts,
        rotor_transitions
    )

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

    size = sum(stator_int[x,y] for x,y in stator_cells)
    model.Minimize(size)

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

        # Only print the pattern if:
        # * it has a smaller population than the original, or
        # * the original stator contains a cell outside the search area, or
        # * there is a cell in the original stator that is forced off in the solution, or
        # * there is a cell not in the original stator that is forced on in the solution.
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
        verbose_print("UNKNOWN: no solutions were found, but the model was not proven infeasible.")
    elif status == cp_model.INFEASIBLE:
        verbose_print("INFEASIBLE: no solution possible.")
    verbose_print()


if __name__ == "__main__":
    main(parse_arguments())
