import argparse
import re


def non_negative_int(value):
    i = int(value)
    if i < 0:
        raise argparse.ArgumentTypeError(f"Value must be 0 or greater. Got: {value}")
    return i


def parse_arguments():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="A program to optimize the stator of patterns "
                    "in Life-like cellular automata.",
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
        type=non_negative_int,
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
        "-s", "--symmetry",
        type=str.upper,
        choices=["C1", "C2", "C4", "D2-", "D2|",
                 "D2/", "D2\\", "D4+", "D4x", "D8"],
        default="C1",
        metavar="SYMMETRY",
        help='Force the stator to have the given symmetry type.\n'
             'Available types are:\n'
             '  "C1" (default), "C2", "C4", "D2-", "D2|",\n'
             '  "D2/", "D2\\", "D4+", "D4x", and "D8".\n'
             'Symmetry is applied relative to the center of the\n'
             'search area. Always enclose the symmetry type in\n'
             'double quotes (""). For D2| and D2\\ symmetries you\n'
             'may need to type "D2\\|" and "D2\\\\" respectively.'
    )
    parser.add_argument(
        "-d", "--distance",
        type=non_negative_int,
        default=None,
        help="Force stator cells to be within this distance of the\n"
             "rotor using the taxicab metric. For still life searches\n"
             "the distance is applied to the forced cells, rather\n"
             "than the rotor."
    )
    parser.add_argument(
        "-b", "--boundary",
        type=str.lower,
        choices=["off", "any"],
        default="off",
        help="The state of the cells at the boundary of the search\n"
             "area (default: 'off'). If 'any' is chosen, then the CA\n"
             "rules will not be applied at the boundary."
    )
    parser.add_argument(
        "-o", "--objectives",
        type=str.lower,
        nargs="*",
        choices=["min_pop", "max_pop", "area", "width", "height", "symmetry",
                 "left", "right", "top", "bottom", "nw", "ne", "sw", "se",
                 "diag", "back_diag"],
        default=["min_pop"],
        metavar="OBJECTIVES",
        help="A list of properties to optimize ordered by priority.\n"
             "Available options are:\n"
             "  min_pop, max_pop, area, width, height, symmetry,\n"
             "  left, right, top, bottom, nw, ne, sw, se, diag,\n"
             "  and back_diag."
    )
    parser.add_argument(
        "--solution_only",
        action="store_true",
        help="Only print the solution to the optimization problem.\n"
             "If there is no solution or if the input pattern is\n"
             "already optimal, print nothing."
    )
    return parser.parse_args()


# Returns the rulestring and rotor grid if the input string can
# be read as a rotor descriptor. Otherwise returns None, None.
def parse_rotor_descriptor(input_string):
    RE_ROTOR = re.compile(r"""
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
    re_match = RE_ROTOR.match(input_string)
    if not re_match:
        return None, None
    rulestring, height, width, grid = re_match.groups()
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


def read_file(file_name):
    try:
        with open(file_name, 'r') as file:
            rle = file.read()
    except FileNotFoundError:
        raise FileNotFoundError("Pattern file not found.")

    return rle
