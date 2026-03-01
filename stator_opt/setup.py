import os
import tempfile

from .rules import get_rulestring


def create_ruletable_session(rulestring, ruletable, lifelib):
    with tempfile.NamedTemporaryFile(   mode='w',
                                        delete=False,
                                        suffix='.txt'
                                    ) as ruletable_file:
            ruletable_file.write(ruletable)
            ruletable_path = ruletable_file.name
    try:
        sess = lifelib.load_rules(
            rulestring,
            "bs8",
            "b12345678s012345678",
            "b1e2-cn3-c4-c5678s012345678",
            ruletable_path
        )
    finally:
        if os.path.exists(ruletable_path):
            os.remove(ruletable_path)

    return sess


def get_pattern(rle, lt, lt5):
    rulestring = get_rulestring(rle, lt)

    # Convert the input file to a history RLE. It's faster to change the
    # cell states by modifying the RLE than to change them in the lifetree.
    initial_pattern = lt5.pattern("", rulestring + "History")
    initial_pattern += lt5.pattern(rle)
    rle = initial_pattern.rle_string()

    # Identify cells to be forced on or off in the solution
    rle_forced_on = rle.translate(str.maketrans("ABC" + "DE" + "F",
                                                "BBB" + "AA" + "."
                                               ))
    rle_forced_off = rle.translate(str.maketrans("ADE" + "BC" + "F",
                                                 "BBB" + "AA" + "."
                                                ))
    forced_on = lt.pattern("", rulestring) + lt5.pattern(rle_forced_on)
    forced_off = lt.pattern("", rulestring) + lt5.pattern(rle_forced_off)

    # Convert all on cells to state 5 and all off cells to state 2.
    # In LifeHistory, State 5 cells will only remain state 5 as long as
    # they never die, so they can be used to detect the initial stator.
    rle = rle.translate(str.maketrans("BD" + "ACE" + "F",
                                      "BB" + "EEE" + "."
                                     ))
    initial_pattern = lt5.pattern(rle)
    if initial_pattern.empty():
        raise ValueError("Input pattern is empty.")

    return initial_pattern, forced_on, forced_off
