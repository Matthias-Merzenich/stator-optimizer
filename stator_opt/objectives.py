from .optimizer import apply_conditional_transform
from .symmetry import has_symmetric_stator

# The symmetry weights were chosen so that
#
# 1. The weight of each symmetry is greater than that of its proper
#    subsymmetries.
# 2. The sum of all subsymmetry weights for each symmetry gives a
#    hierarchy of symmetry preference.
# 3. C2 is favored over bilateral symmetry across a single axis, but
#    is disfavored over bilateral symmetry across two axes.
# 4. C4 is favored over bilateral symmetry across two orthogonal axes.
#
# Conditions 3 and 4 were chosen based on a community vote showing
# that rotational symmetry was more popular than bilateral symmetry.
SYMMETRY_WEIGHTS = {
    "C1": 0,    "C2": 2,    "C4": 6,    "D2-": 1,   "D2|": 1,
    "D2/": 1,   "D2\\": 1,  "D4+": 3,   "D4X": 3,   "D8": 7
}

SYMMETRY_LIST = list(SYMMETRY_WEIGHTS)

OBJECTIVE_ALIAS = {
    "min_pop": "Stator population",
    "max_pop": "Stator population",
    "min_change": "Changed",
    "max_change": "Changed",
    "boundary_pop": "Boundary population",
    "valid_boundary": "Invalid boundary cells",
    "symmetry": "Symmetry",
    "left": "x_min",
    "right": "x_max",
    "top": "y_min",
    "bottom": "y_max",
    "width": "Width",
    "height": "Height",
    "area": "Area",
    "ne": "ne_max",
    "sw": "sw_min",
    "se": "se_max",
    "nw": "nw_min",
    "diag": "Diagonal",
    "back_diag": "Backwards diagonal",
    "left_with_rotor": "x_min (with rotor)",
    "right_with_rotor": "x_max (with rotor)",
    "top_with_rotor": "y_min (with rotor)",
    "bottom_with_rotor": "y_max (with rotor)",
    "width_with_rotor": "Width (with rotor)",
    "height_with_rotor": "Height (with rotor)",
    "area_with_rotor": "Area (with rotor)",
    "ne_with_rotor": "ne_max (with rotor)",
    "sw_with_rotor": "sw_min (with rotor)",
    "se_with_rotor": "se_max (with rotor)",
    "nw_with_rotor": "nw_min (with rotor)",
    "diag_with_rotor": "Diagonal (with rotor)",
    "back_diag_with_rotor": "Backwards diagonal (with rotor)"
}


class Extrema:
    def __init__(self):
        self.x_max = None
        self.x_min = None
        self.y_max = None
        self.y_min = None
        self.width = None
        self.height = None
        self.area = None
        self.ne = None
        self.sw = None
        self.se = None
        self.nw = None
        self.diag = None
        self.back_diag = None


class ObjectiveStats:
    def __init__(self, the_pattern, any_boundary):
        self._status = "Initial"
        if the_pattern.is_history:
            init_symmetries = [symm for symm in SYMMETRY_LIST
                               if has_symmetric_stator(the_pattern, symm)]
            self.stats = {
                "min_pop": the_pattern.initial_stator_on.population,
                "max_pop": the_pattern.initial_stator_on.population,
                "min_change": 0,
                "max_change": 0,
                "symmetry": sum(SYMMETRY_WEIGHTS[symm]
                                for symm in init_symmetries),
                "symmetry_name": (
                    f'"{max(init_symmetries, key=SYMMETRY_WEIGHTS.get)}"'
                ),
                "boundary_pop": (
                    the_pattern.initial_stator_on
                    & the_pattern.stator
                    & any_boundary
                ).population,
                "valid_boundary": 0
            }
            for with_rotor in [True, False]:
                self.stats.update(get_box_stats(the_pattern, with_rotor))
        else:
            self.stats = {
                "min_pop": float("inf"),
                "max_pop": - float("inf"),
                "min_change": float("inf"),
                "max_change": - float("inf"),
                "symmetry": - float("inf"),
                "symmetry_name": '"C1"',
                "boundary_pop": float("inf"),
                "valid_boundary": float("inf"),
                "left": - float("inf"),    "left_with_rotor": - float("inf"),
                "right": float("inf"),     "right_with_rotor": float("inf"),
                "top": - float("inf"),     "top_with_rotor": - float("inf"),
                "bottom": float("inf"),    "bottom_with_rotor": float("inf"),
                "width": float("inf"),     "width_with_rotor": float("inf"),
                "height": float("inf"),    "height_with_rotor": float("inf"),
                "area": float("inf"),      "area_with_rotor": float("inf"),
                "ne": float("inf"),        "ne_with_rotor": float("inf"),
                "sw": - float("inf"),      "sw_with_rotor": - float("inf"),
                "se": float("inf"),        "se_with_rotor": float("inf"),
                "nw": - float("inf"),      "nw_with_rotor": - float("inf"),
                "diag": float("inf"),      "diag_with_rotor": float("inf"),
                "back_diag": float("inf"), "back_diag_with_rotor": float("inf")
            }

    def store_final_stats(
        self, solver, stator_vars, box_vars, symm_vars,
        change_var, boundary_pop_var, invalid_boundary_sum
    ):
        self._status = "Final"
        final_symmetries = [
            symm for symm in SYMMETRY_LIST if solver.Value(symm_vars[symm])
        ]
        if not final_symmetries:
            final_symmetries = ["C1"]
        final_stator_pop = sum(solver.Value(v) for v in stator_vars.values())
        self.stats = {
            "min_pop": final_stator_pop,
            "max_pop": final_stator_pop,
            "min_change": solver.Value(change_var),
            "max_change": solver.Value(change_var),
            "symmetry": sum(SYMMETRY_WEIGHTS[symm]
                            for symm in final_symmetries),
            "symmetry_name": (
                f'"{max(final_symmetries, key=SYMMETRY_WEIGHTS.get)}"'
            ),
            "boundary_pop": solver.Value(boundary_pop_var),
            "valid_boundary": solver.Value(invalid_boundary_sum)
        }
        for key in box_vars:
            self.stats[key] = solver.Value(box_vars[key])

    def get_stats_string(self, objectives):
        lines = [f"\n{self._status} objective values:"]
        for obj in objectives:
            if obj == "symmetry":
                lines.append(f"  Symmetry: {self.stats['symmetry_name']}")
            else:
                lines.append(f"  {OBJECTIVE_ALIAS[obj]}: {self.stats[obj]}")
        return "\n".join(lines) + "\n"


# Given a pattern, return an Extrema object containing various extreme values.
def get_extrema(patt):
    cell_array = patt.coords()
    extrema = Extrema()
    extrema.x_max, extrema.y_max = cell_array.max(axis=0).tolist()
    extrema.x_min, extrema.y_min = cell_array.min(axis=0).tolist()
    cell_array_diff = cell_array[:, 0] - cell_array[:, 1]
    extrema.ne = cell_array_diff.max().tolist()
    extrema.sw = cell_array_diff.min().tolist()
    extrema.se = cell_array.sum(axis=1).max().tolist()
    extrema.nw = cell_array.sum(axis=1).min().tolist()
    return extrema


def get_box_stats(the_pattern, with_rotor):
    if with_rotor:
        suffix = "_with_rotor"
        max_pattern = the_pattern.stator + the_pattern.clipped_rotor
        init_pattern = the_pattern.initial_stator_on + the_pattern.clipped_rotor
    else:
        suffix = ""
        max_pattern = the_pattern.stator
        init_pattern = the_pattern.initial_stator_on

    if init_pattern.population == 0:
        init_pattern = max_pattern

    init = get_extrema(init_pattern)

    init.width = init.x_max - init.x_min + 1
    init.height = init.y_max - init.y_min + 1
    init.diag = init.ne - init.sw + 1
    init.back_diag = init.se - init.nw + 1

    stats = {
        "left" + suffix: init.x_min,
        "right" + suffix: init.x_max,
        "top" + suffix: init.y_min,
        "bottom" + suffix: init.y_max,
        "width" + suffix: init.width,
        "height" + suffix: init.height,
        "area" + suffix: init.width * init.height,
        "ne" + suffix: init.ne,
        "sw" + suffix: init.sw,
        "se" + suffix: init.se,
        "nw" + suffix: init.nw,
        "diag" + suffix: init.diag,
        "back_diag" + suffix: init.back_diag
    }

    return stats


def apply_box_objectives(
    model, the_pattern, objectives, stator_vars
):
    # The construction of variables for objectives with rotor is essentially
    # the same as for objectives without rotor. Order is important, because
    # we later use some of the values calculated when with_rotor is False.
    box_vars = {}
    objective_dict = {}
    for with_rotor in [True, False]:
        if with_rotor:
            wr = "_with_rotor"  # Suffix for objective names
            max_pattern = the_pattern.base_stator + the_pattern.clipped_rotor
        else:
            wr = ""
            max_pattern = the_pattern.base_stator

        extreme = get_extrema(max_pattern)

        extreme.width = extreme.x_max - extreme.x_min + 1
        extreme.height = extreme.y_max - extreme.y_min + 1
        extreme.area = extreme.width * extreme.height
        extreme.diag = extreme.ne - extreme.sw + 1
        extreme.back_diag = extreme.se - extreme.nw + 1

        x_max_var = model.NewIntVar(extreme.x_min, extreme.x_max, "x_max" + wr)
        x_min_var = model.NewIntVar(extreme.x_min, extreme.x_max, "x_min" + wr)
        y_max_var = model.NewIntVar(extreme.y_min, extreme.y_max, "y_max" + wr)
        y_min_var = model.NewIntVar(extreme.y_min, extreme.y_max, "y_min" + wr)
        width_var = model.NewIntVar(0, extreme.width, "width" + wr)
        height_var = model.NewIntVar(0, extreme.height, "height" + wr)
        area_var = model.NewIntVar(0, extreme.area, "area" + wr)
        ne_var = model.NewIntVar(extreme.sw, extreme.ne, "ne" + wr)
        sw_var = model.NewIntVar(extreme.sw, extreme.ne, "sw" + wr)
        se_var = model.NewIntVar(extreme.nw, extreme.se, "se" + wr)
        nw_var = model.NewIntVar(extreme.nw, extreme.se, "nw" + wr)
        diag_var = model.NewIntVar(0, extreme.diag, "diag" + wr)
        back_diag_var = model.NewIntVar(0, extreme.back_diag, "back_diag" + wr)

        box_vars.update({
            "left" + wr: x_min_var,
            "right" + wr: x_max_var,
            "top" + wr: y_min_var,
            "bottom" + wr: y_max_var,
            "width" + wr: width_var,
            "height" + wr: height_var,
            "area" + wr: area_var,
            "ne" + wr: ne_var,
            "sw" + wr: sw_var,
            "se" + wr: se_var,
            "nw" + wr: nw_var,
            "diag" + wr: diag_var,
            "back_diag" + wr: back_diag_var
        })

        objective_dict.update({
            "left"+wr:     [box_vars["left"+wr],     "max", extreme.width],
            "right"+wr:    [box_vars["right"+wr],    "min", extreme.width],
            "top"+wr:      [box_vars["top"+wr],      "max", extreme.height],
            "bottom"+wr:   [box_vars["bottom"+wr],   "min", extreme.height],
            "width"+wr:    [box_vars["width"+wr],    "min", extreme.width],
            "height"+wr:   [box_vars["height"+wr],   "min", extreme.height],
            "area"+wr:     [box_vars["area"+wr],     "min", extreme.area],
            "ne"+wr:       [box_vars["ne"+wr],       "min", extreme.diag],
            "sw"+wr:       [box_vars["sw"+wr],       "max", extreme.diag],
            "se"+wr:       [box_vars["se"+wr],       "min", extreme.back_diag],
            "nw"+wr:       [box_vars["nw"+wr],       "max", extreme.back_diag],
            "diag"+wr:     [box_vars["diag"+wr],     "min", extreme.diag],
            "back_diag"+wr:[box_vars["back_diag"+wr],"min", extreme.back_diag]
        })

    objective_set = set(objectives)
    constraint_needed = {
        "right": {"right", "width", "area", "right_with_rotor",
                  "width_with_rotor", "area_with_rotor"} & objective_set,
        "left": {"left", "width", "area", "left_with_rotor",
                 "width_with_rotor", "area_with_rotor"} & objective_set,
        "bottom": {"bottom", "height", "area", "bottom_with_rotor",
                   "height_with_rotor", "area_with_rotor"} & objective_set,
        "top": {"top", "height", "area", "top_with_rotor",
                "height_with_rotor", "area_with_rotor"} & objective_set,
        "ne": {"ne", "diag", "ne_with_rotor",
               "diag_with_rotor"} & objective_set,
        "sw": {"sw", "diag", "sw_with_rotor",
               "diag_with_rotor"} & objective_set,
        "se": {"se", "back_diag", "se_with_rotor",
               "back_diag_with_rotor"} & objective_set,
        "nw": {"nw", "back_diag", "nw_with_rotor",
               "back_diag_with_rotor"} & objective_set
    }

    def cell_set(func, extremum):
        return [
            func(x,y) * stator_vars[x,y] + (1 - stator_vars[x,y]) * extremum
            for x,y in stator_vars
        ]

    # Add constraints to the model so that the variables actually do what we
    # want them to. We can use our named tuple (extreme) and direct variable
    # names here, because we just used with_rotor=False in the above loop,
    # meaning the variables and extremes are for box variables without rotor.
    if constraint_needed["right"]:
        model.AddMaxEquality(x_max_var, cell_set(lambda x, y: x, extreme.x_min))
    if constraint_needed["left"]:
        model.AddMinEquality(x_min_var, cell_set(lambda x, y: x, extreme.x_max))
    if constraint_needed["bottom"]:
        model.AddMaxEquality(y_max_var, cell_set(lambda x, y: y, extreme.y_min))
    if constraint_needed["top"]:
        model.AddMinEquality(y_min_var, cell_set(lambda x, y: y, extreme.y_max))
    if constraint_needed["ne"]:
        model.AddMaxEquality(ne_var, cell_set(lambda x, y: x - y, extreme.sw))
    if constraint_needed["sw"]:
        model.AddMinEquality(sw_var, cell_set(lambda x, y: x - y, extreme.ne))
    if constraint_needed["se"]:
        model.AddMaxEquality(se_var, cell_set(lambda x, y: x + y, extreme.nw))
    if constraint_needed["nw"]:
        model.AddMinEquality(nw_var, cell_set(lambda x, y: x + y, extreme.se))
    model.Add(width_var == x_max_var - x_min_var + 1)
    model.Add(height_var == y_max_var - y_min_var + 1)
    if {"area"} & objective_set:
        model.AddMultiplicationEquality(area_var, [width_var, height_var])
    model.Add(diag_var == ne_var - sw_var + 1)
    model.Add(back_diag_var == se_var - nw_var + 1)

    if the_pattern.clipped_rotor.nonempty():
        rotor_extreme = get_extrema(the_pattern.clipped_rotor)
        model.AddMaxEquality(box_vars["right_with_rotor"],
                             [x_max_var, rotor_extreme.x_max])
        model.AddMinEquality(box_vars["left_with_rotor"],
                             [x_min_var, rotor_extreme.x_min])
        model.AddMaxEquality(box_vars["bottom_with_rotor"],
                             [y_max_var, rotor_extreme.y_max])
        model.AddMinEquality(box_vars["top_with_rotor"],
                             [y_min_var, rotor_extreme.y_min])
        model.AddMaxEquality(box_vars["ne_with_rotor"],
                             [ne_var, rotor_extreme.ne])
        model.AddMinEquality(box_vars["sw_with_rotor"],
                             [sw_var, rotor_extreme.sw])
        model.AddMaxEquality(box_vars["se_with_rotor"],
                             [se_var, rotor_extreme.se])
        model.AddMinEquality(box_vars["nw_with_rotor"],
                             [nw_var, rotor_extreme.nw])
        model.Add(
            box_vars["width_with_rotor"] == (
                box_vars["right_with_rotor"] - box_vars["left_with_rotor"] + 1
            )
        )
        model.Add(
            box_vars["height_with_rotor"] == (
                box_vars["bottom_with_rotor"] - box_vars["top_with_rotor"] + 1
            )
        )
        if {"area_with_rotor"} & objective_set:
            model.AddMultiplicationEquality(area_var,
                [box_vars["width_with_rotor"], box_vars["height_with_rotor"]]
            )
        model.Add(
            box_vars["diag_with_rotor"] == (
                box_vars["ne_with_rotor"] - box_vars["sw_with_rotor"] + 1
            )
        )
        model.Add(
            box_vars["back_diag_with_rotor"] == (
                box_vars["se_with_rotor"] - box_vars["nw_with_rotor"] + 1
            )
        )

    return box_vars, objective_dict


def apply_symmetry_objective(
    model, objectives, symmetry, stator_vars, stator_array
):
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

    for symm in SYMMETRY_LIST:
        symm_bool[symm] = model.NewBoolVar(f"symm_bool_{symm}")
        for trans in FORCED_TRANSFORMATIONS[symm]:
            model.AddImplication(symm_bool[symm], trans_bool[trans])

    transformations = list(FORCED_TRANSFORMATIONS[symmetry])
    if "symmetry" in objectives:
        transformations = TRANSFORMATION_LIST

    for transformation in transformations:
        apply_conditional_transform(
            model,
            transformation,
            trans_bool,
            stator_vars,
            stator_array
        )

    # Cause transformations to be applied unconditionally
    for trans in FORCED_TRANSFORMATIONS[symmetry]:
        model.Add(trans_bool[trans] == 1)

    return symm_bool


def get_objective_expression(
    the_pattern, objective_stats, objective_dict, objectives
):
    current_weight = 1
    total_objective = 0
    init_objective = 0
    objective_upper_bound = 0

    for key in reversed(objectives):
        var, objective_type, max_reduction_amount = objective_dict[key]
        min_or_max = 1 if objective_type == "min" else -1
        init_value = objective_stats.stats[key]
        total_objective += var * min_or_max * current_weight
        init_objective += init_value * min_or_max * current_weight
        objective_upper_bound += max_reduction_amount * current_weight
        current_weight *= max_reduction_amount + 1

    if not the_pattern.is_history:
        init_objective = float("inf")

    if abs(objective_upper_bound) > 2**63 - 1:
        raise OverflowError("Objective expression exceeds 64-bit limit.")

    return total_objective, init_objective
