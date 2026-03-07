import numpy as np


# Treat coordinate tuples as single objects
def coordinate_view(arr):
    arr = np.ascontiguousarray(arr)
    return arr.view(np.dtype((np.void, arr.dtype.itemsize * arr.shape[1])))


def coord_intersect(coords_1, coords_2):
    intersection = np.intersect1d(
        coordinate_view(coords_1), coordinate_view(coords_2)
    )
    return intersection.view(coords_1.dtype).reshape(-1, coords_1.shape[1])


def transform(coords, mode, center):
    # coords[:, [1, 0]] flips coordinates across y=x
    # coords * [a, b] multiplies componentwise
    transform_dict = {
        "flip_x": lambda c: c * [-1, 1],
        "flip_y": lambda c: c * [1, -1],
        "flip_diag": lambda c: c[:, [1, 0]],
        "flip_reverse_diag": lambda c: c[:, [1, 0]] * [-1, -1],
        "rotate_90": lambda c: c[:, [1, 0]] * [-1, 1],
        "rotate_180": lambda c: c * [-1, -1],
        "rotate_270": lambda c: c[:, [1, 0]] * [1, -1]
    }
    return (transform_dict[mode](coords*2 - center) + center) // 2


def symmetrize_with_center(coords, symmetry, center):
    if (
        symmetry in ["C4", "D2/", "D2\\", "D4X", "D8"]
        and (center[0] % 2 != center[1] % 2)
    ):
        return np.empty((0, 2), dtype=np.int64), center

    new_coords = coords.copy()
    if symmetry in [None, "C1"]:
        return new_coords, center

    TRANSFORMS = {
        "C1": [],
        "C2": ["rotate_180"],
        "C4": ["rotate_90", "rotate_180"],
        "D2-": ["flip_y"],
        "D2|": ["flip_x"],
        "D2/": ["flip_diag"],
        "D2\\": ["flip_reverse_diag"],
        "D4+": ["flip_x", "flip_y"],
        "D4X": ["flip_diag", "flip_reverse_diag"],
        "D8": ["flip_x", "flip_y", "flip_diag"]
    }
    for mode in TRANSFORMS[symmetry]:
        new_coords = coord_intersect(
            new_coords, transform(new_coords, mode, center)
        )
    return new_coords, center


def symmetrize(coords, symmetry):
    # This is double the actual center of the coordinates.
    # The purpose of this is to avoid half-integer centers.
    center = coords.min(axis=0) + coords.max(axis=0)
    return symmetrize_with_center(coords, symmetry, center)


# This checks if the stator is symmetric about the center
# of the search area, not the center of the original pattern.
def has_symmetric_stator(pattern, symmetry):
    stator_array = pattern.stator.coords()
    center = stator_array.min(axis=0) + stator_array.max(axis=0)
    on_stator_array = pattern.initial_stator_on.coords()
    symmetric_stator_array, _ = symmetrize_with_center(
        on_stator_array, symmetry, center
    )
    on_stator_cells = set(map(tuple, on_stator_array.tolist()))
    symmetric_stator_cells = set(map(tuple, symmetric_stator_array.tolist()))
    return (on_stator_cells == symmetric_stator_cells)
