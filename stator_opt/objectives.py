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
    "symmetry": "Symmetry"
}

class ObjectiveStats:
    def __init__(self, the_pattern):
        self._status = "Initial"
        stator_array = the_pattern.stator.coords()
        max_x, max_y = stator_array.max(axis=0).tolist()
        min_x, min_y = stator_array.min(axis=0).tolist()
        stator_array_diff = stator_array[:, 0] - stator_array[:, 1]
        max_ne = stator_array_diff.max().tolist()
        min_sw = stator_array_diff.min().tolist()
        max_se = stator_array.sum(axis=1).max().tolist()
        min_nw = stator_array.sum(axis=1).min().tolist()

        if the_pattern.initial.getrule().endswith("History"):
            init_stator_array = the_pattern.initial_stator_on.coords()
            if init_stator_array.size == 0:
                init_x_max = max_x
                init_x_min = min_x
                init_y_max = max_y
                init_y_min = min_y
                init_ne = max_ne
                init_sw = min_sw
                init_se = max_se
                init_nw = min_nw
            else:
                init_x_max, init_y_max = init_stator_array.max(axis=0).tolist()
                init_x_min, init_y_min = init_stator_array.min(axis=0).tolist()
                init_stator_array_diff = (
                    init_stator_array[:, 0] - init_stator_array[:, 1]
                )
                init_ne = init_stator_array_diff.max().tolist()
                init_sw = init_stator_array_diff.min().tolist()
                init_se = init_stator_array.sum(axis=1).max().tolist()
                init_nw = init_stator_array.sum(axis=1).min().tolist()
            init_diag = init_ne - init_sw + 1
            init_back_diag = init_se - init_nw + 1
            init_width = init_x_max - init_x_min + 1
            init_height = init_y_max - init_y_min + 1
            init_symmetries = [symm for symm in SYMMETRY_LIST
                               if has_symmetric_stator(the_pattern, symm)]
            self.stats = {
                "min_pop": the_pattern.initial_stator_on.population,
                "max_pop": the_pattern.initial_stator_on.population,
                "left": init_x_min,
                "right": init_x_max,
                "top": init_y_min,
                "bottom": init_y_max,
                "width": init_width,
                "height": init_height,
                "area": init_width * init_height,
                "ne": init_ne,
                "sw": init_sw,
                "se": init_se,
                "nw": init_nw,
                "diag": init_diag,
                "back_diag": init_back_diag,
                "symmetry": sum(SYMMETRY_WEIGHTS[symm]
                                for symm in init_symmetries),
                "symmetry_name": (
                    f'"{max(init_symmetries, key=SYMMETRY_WEIGHTS.get)}"'
                )
            }
        else:
            self.stats = {
                "min_pop": float("inf"),
                "max_pop": - float("inf"),
                "left": - float("inf"),
                "right": float("inf"),
                "top": - float("inf"),
                "bottom": float("inf"),
                "width": float("inf"),
                "height": float("inf"),
                "area": float("inf"),
                "ne": float("inf"),
                "sw": - float("inf"),
                "se": float("inf"),
                "nw": - float("inf"),
                "diag": float("inf"),
                "back_diag": float("inf"),
                "symmetry": - float("inf"),
                "symmetry_name": '"C1"'
            }

    def store_final_stats(self, solver, stator_vars, box_vars, symm_vars):
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
            "symmetry": sum(SYMMETRY_WEIGHTS[symm]
                            for symm in final_symmetries),
            "symmetry_name": (
                f'"{max(final_symmetries, key=SYMMETRY_WEIGHTS.get)}"'
            )
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


def apply_box_objectives(
    model, the_pattern, objectives, stator_array, stator_bool
):
    max_x, max_y = stator_array.max(axis=0).tolist()
    min_x, min_y = stator_array.min(axis=0).tolist()
    max_width = max_x - min_x + 1
    max_height = max_y - min_y + 1
    max_area = max_width * max_height
    stator_array_diff = stator_array[:, 0] - stator_array[:, 1]
    max_ne = stator_array_diff.max().tolist()
    min_sw = stator_array_diff.min().tolist()
    max_se = stator_array.sum(axis=1).max().tolist()
    min_nw = stator_array.sum(axis=1).min().tolist()
    max_diag = max_ne - min_sw + 1
    max_back_diag = max_se - min_nw + 1

    x_max_var = model.NewIntVar(min_x, max_x, "x_max")
    x_min_var = model.NewIntVar(min_x, max_x, "x_min")
    y_max_var = model.NewIntVar(min_y, max_y, "y_max")
    y_min_var = model.NewIntVar(min_y, max_y, "y_min")
    width_var = model.NewIntVar(0, max_width, "width")
    height_var = model.NewIntVar(0, max_height, "height")
    area_var = model.NewIntVar(0, max_area, "area")
    ne_var = model.NewIntVar(min_sw, max_ne, "ne")
    sw_var = model.NewIntVar(min_sw, max_ne, "sw")
    se_var = model.NewIntVar(min_nw, max_se, "se")
    nw_var = model.NewIntVar(min_nw, max_se, "nw")
    diag_var = model.NewIntVar(0, max_diag, "diag")
    back_diag_var = model.NewIntVar(0, max_back_diag, "back_diag")

    box_vars = {
        "left": x_min_var,
        "right": x_max_var,
        "top": y_min_var,
        "bottom": y_max_var,
        "width": width_var,
        "height": height_var,
        "area": area_var,
        "ne": ne_var,
        "sw": sw_var,
        "se": se_var,
        "nw": nw_var,
        "diag": diag_var,
        "back_diag": back_diag_var
    }

    objective_set = set(objectives)
    constraint_needed = {
        "right": {"right", "width", "area"} & objective_set,
        "left": {"left", "width", "area"} & objective_set,
        "bottom": {"bottom", "height", "area"} & objective_set,
        "top": {"top", "height", "area"} & objective_set,
        "ne": {"ne", "diag"} & objective_set,
        "sw": {"sw", "diag"} & objective_set,
        "se": {"se", "back_diag"} & objective_set,
        "nw": {"nw", "back_diag"} & objective_set
    }

    def cell_set(func, extremum):
        return [
            func(x,y) * stator_bool[x,y] + (1 - stator_bool[x,y]) * extremum
            for x,y in stator_bool
        ]

    if constraint_needed["right"]:
        model.AddMaxEquality(x_max_var, cell_set(lambda x, y: x, min_x))
    if constraint_needed["left"]:
        model.AddMinEquality(x_min_var, cell_set(lambda x, y: x, max_x))
    if constraint_needed["bottom"]:
        model.AddMaxEquality(y_max_var, cell_set(lambda x, y: y, min_y))
    if constraint_needed["top"]:
        model.AddMinEquality(y_min_var, cell_set(lambda x, y: y, max_y))
    if constraint_needed["ne"]:
        model.AddMaxEquality(ne_var, cell_set(lambda x, y: x - y, min_sw))
    if constraint_needed["sw"]:
        model.AddMinEquality(sw_var, cell_set(lambda x, y: x - y, max_ne))
    if constraint_needed["se"]:
        model.AddMaxEquality(se_var, cell_set(lambda x, y: x + y, min_nw))
    if constraint_needed["nw"]:
        model.AddMinEquality(nw_var, cell_set(lambda x, y: x + y, max_se))
    model.Add(width_var == x_max_var - x_min_var + 1)
    model.Add(height_var == y_max_var - y_min_var + 1)
    if {"area"} & objective_set:
        model.AddMultiplicationEquality(area_var, [width_var, height_var])
    model.Add(diag_var == ne_var - sw_var + 1)
    model.Add(back_diag_var == se_var - nw_var + 1)

    objective_dict = {
        "left":      [box_vars["left"],      "max",     max_width],
        "right":     [box_vars["right"],     "min",     max_width],
        "top":       [box_vars["top"],       "max",     max_height],
        "bottom":    [box_vars["bottom"],    "min",     max_height],
        "width":     [box_vars["width"],     "min",     max_width],
        "height":    [box_vars["height"],    "min",     max_height],
        "area":      [box_vars["area"],      "min",     max_area],
        "ne":        [box_vars["ne"],        "min",     max_diag],
        "sw":        [box_vars["sw"],        "max",     max_diag],
        "se":        [box_vars["se"],        "min",     max_back_diag],
        "nw":        [box_vars["nw"],        "max",     max_back_diag],
        "diag":      [box_vars["diag"],      "min",     max_diag],
        "back_diag": [box_vars["back_diag"], "min",     max_back_diag]
    }

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

    if not the_pattern.initial.getrule().endswith("History"):
        init_objective = float("inf")

    if abs(objective_upper_bound) > 2**63 - 1:
        raise OverflowError("Objective expression exceeds 64-bit limit.")

    return total_objective, init_objective
