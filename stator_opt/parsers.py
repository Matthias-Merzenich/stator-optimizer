import argparse
import re
import textwrap


# Custom argparse formatter class that removes redundancy
# and gives more control over help text.
class CustomFormatter(argparse.HelpFormatter):
    # Wrap help text lines while still respecting manual newlines.
    def _split_lines(self, text, width):
        MAX_WIDTH = 90
        effective_width = min(width, MAX_WIDTH)
        lines = []
        for paragraph in text.splitlines():
            if paragraph.strip() == "":
                # Preserve blank lines.
                lines.append("")
            else:
                # Wrap paragraph if it's too long.
                lines.extend(textwrap.wrap(paragraph,
                                           width=effective_width,
                                           replace_whitespace=False))
        return lines

    # Print metavars as written, without additional ellipses.
    def _format_args(self, action, default_metavar):
        if action.metavar is not None:
            if isinstance(action.metavar, tuple):
                return " ".join(action.metavar)
            return action.metavar
        return super()._format_args(action, default_metavar)

    # List all option strings before the metavar to reduce
    # redundancy and keep the help lines short.
    def _format_action_invocation(self, action):
        # Use default for positional arguments
        if not action.option_strings:
            return super()._format_action_invocation(action)

        if action.nargs != 0:
            args_string = self._format_args(action, action.dest.upper())
            options = ", ".join(action.option_strings)
            return f"{options} {args_string}"

        return ", ".join(action.option_strings)


# Adjustments should have the form (int, int, int, int, [{"off","any"}])
class AddAdjustments(argparse.Action):
    def __call__(self, parser, namespace, items, option_string=None):
        if not (4 <= len(items) <= 5):
            raise argparse.ArgumentTypeError(
                f"option {option_string}: "
                "Argument must have 4 or 5 elements."
            )
        try:
            adjustments = [int(items[i]) for i in range(4)]
        except ValueError:
            raise argparse.ArgumentTypeError(
                f"option {option_string}: "
                f"First four elements must be integers. "
                f"Got: '{items[:4]}'."
            )
        boundary = "off"
        if len(items) == 5:
            boundary = items[4].lower()
            if boundary not in ["off", "any"]:
                raise argparse.ArgumentTypeError(
                    f"option {option_string}: "
                    f"Fifth element must be 'off' or 'any'. "
                    f"Got: '{items[4]}'."
                )

        # Update values
        current_options = getattr(namespace, self.dest)
        current_options["adjustments"] = adjustments
        current_options["boundary"] = boundary
        setattr(namespace, self.dest, current_options)


# Distances should have the form ({"rotor","stator"}, int>=0, [{"off","any"}])
class AddDistance(argparse.Action):
    def __call__(self, parser, namespace, items, option_string=None):
        if not (2 <= len(items) <= 3):
            raise argparse.ArgumentTypeError(
                f"option {option_string}: "
                "Argument must have 2 or 3 elements."
            )
        source = items[0].lower()
        if source not in ["rotor", "stator"]:
            raise argparse.ArgumentTypeError(
                f"option {option_string}: "
                f"First element must be 'rotor' or 'stator'. "
                f"Got '{items[0]}'."
            )
        try:
            distance = non_negative_int(items[1])
        except argparse.ArgumentTypeError:
            raise argparse.ArgumentTypeError(
                f"option {option_string}: "
                f"Second element must be a nonnegative integer. "
                f"Got: '{items[1]}'."
            )
        boundary = "off"
        if len(items) == 3:
            boundary = items[2].lower()
            if boundary not in ["off", "any"]:
                raise argparse.ArgumentTypeError(
                    f"option {option_string}: "
                    f"Third element must be 'off' or 'any'. "
                    f"Got: '{items[2]}'."
                )

        # Update values
        current_options = getattr(namespace, self.dest)
        current_options[source]["distance"] = distance
        current_options[source]["boundary"] = boundary
        setattr(namespace, self.dest, current_options)


def non_negative_int(value):
    try:
        i = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"Value must be an integer. Got: {value}")
    if i < 0:
        raise argparse.ArgumentTypeError(
            f"Value must be 0 or greater. Got: {value}"
        )
    return i


def parse_arguments():
    parser = argparse.ArgumentParser(
        formatter_class=CustomFormatter,
        description="A program to optimize the stator of patterns "
                    "in Life-like cellular automata.",
        add_help=False
    )
    parser.add_argument(
        "input_file",
        help="File containing the pattern to be optimized. When the "
             "pattern is saved using a history rule, you can specify "
             "stator cell states in the solution as follows:\n"
             "  ------------------------------------------------\n"
             "   history state | input state | solution state\n"
             "  ---------------+-------------+------------------\n"
             "   0 (black)     | off         | either on or off\n"
             "   1 (green)     | on          | either on or off\n"
             "   2 (blue)      | off         | off\n"
             "   3 (white)     | on          | off\n"
             "   4 (red)       | off         | on\n"
             "   5 (yellow)    | on          | on\n"
             "   6 (gray)      | off         | either on or off\n"
             "This only applies to stator cells. If a cell is "
             "determined to be part of the rotor, then its state will "
             "not be changed in the solution.\n"
             "Accepted input formats are RLE (.rle), macrocell (.mc), "
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
        nargs="+",
        action=AddAdjustments,
        metavar=("LEFT", "RIGHT", "TOP", "BOTTOM", "[{off,any}]"),
        help="Expand the search box by the given distances; negative "
             "values contract it. An optional boundary state can also "
             "be specified; if set to 'any', then the CA rules will "
             "not be applied at the boundary."
    )
    parser.add_argument(
        "-d", "--distance",
        nargs="+",
        action=AddDistance,
        metavar=("{rotor,stator}", "DISTANCE", "[{off,any}]"),
        help="Force stator cells to be within the given distance of "
             "either the rotor or the initial stator using the "
             "taxicab metric. To set distances from both the rotor "
             "and the stator, specify the option twice. For still "
             "life searches, the distance is measured from the forced "
             "cells, rather than the rotor. An optional boundary "
             "state can also be specified; if set to 'any', then the "
             "CA rules will not be applied at the boundary."
    )
    parser.add_argument(
        "-s", "--symmetry",
        type=str.upper,
        choices=["C1", "C2", "C4", "D2-", "D2|",
                 "D2/", "D2\\", "D4+", "D4x", "D8"],
        metavar="SYMMETRY",
        help='Force the stator to have the given symmetry type.\n'
             'Available types:\n'
             '  "C1" (default), "C2", "C4", "D2-", "D2|",\n'
             '  "D2/", "D2\\", "D4+", "D4x", and "D8".\n'
             'Symmetry is applied relative to the center of the '
             'search area. Always enclose the symmetry type in double '
             'quotes (""); for D2| and D2\\ symmetries you may need to '
             'type "D2\\|" and "D2\\\\", respectively.'
    )
    parser.add_argument(
        "-o", "--optimize",
        type=str.lower,
        nargs="*",
        choices=["min_pop", "max_pop", "area", "width", "height", "symmetry",
                 "left", "right", "top", "bottom", "nw", "ne", "sw", "se",
                 "diag", "back_diag", "area_with_rotor", "width_with_rotor",
                 "height_with_rotor", "left_with_rotor", "right_with_rotor",
                 "top_with_rotor", "bottom_with_rotor", "nw_with_rotor",
                 "ne_with_rotor", "sw_with_rotor", "se_with_rotor",
                 "diag_with_rotor", "back_diag_with_rotor"],
        metavar="[OBJECTIVES ...]",
        help="A list of properties to optimize ordered by priority.\n"
             "Available properties:\n"
             "  min_pop, max_pop, area, width, height, symmetry,\n"
             "  left, right, top, bottom, nw, ne, sw, se, diag,\n"
             "  and back_diag.\n"
             "The calculations of these properties only use the "
             "stator. Append '_with_rotor' to also use the clipped "
             "rotor from the box specified by the -c option."
    )
    parser.add_argument(
        "-c", "--clip_rotor",
        type=int,
        nargs=4,
        metavar=("LEFT", "RIGHT", "TOP", "BOTTOM"),
        help="Expand or contract the clipped rotor box. The initial "
            "box is the bounding box of the input pattern. This "
            "option is only used when '_with_rotor' objectives are "
            "specified by the -o option."
    )
    parser.add_argument(
        "--solution_only",
        action="store_true",
        help="Only print the solution. If no solution exists or if "
             "the input pattern is already optimal, print nothing."
    )
    parser.set_defaults(
        adjust={"adjustments": [0, 0, 0, 0], "boundary": "off"},
        symmetry="C1",
        distance={
            "rotor": {"distance": None, "boundary": "off"},
            "stator": {"distance": None, "boundary": "off"}
        },
        boundary="off",
        optimize=["min_pop"],
        clip_rotor=[0, 0, 0, 0])
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
        raise FileNotFoundError(f"Pattern file not found: {file_name}")

    return rle
