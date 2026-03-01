import re


def get_rulestring_from_pattern(pattern):
    rulestring = pattern.getrule()
    if rulestring.endswith("History"):
        rulestring = rulestring[:-len("History")]
    return rulestring


def get_rulestring(rle, lt):
    return get_rulestring_from_pattern(lt.pattern(rle))


def get_rule(rulestring):
    match = re.match(r"b([2-8]*)s([0-8]*)$", rulestring)
    if not match:
        raise ValueError(
            f"Invalid rule '{rulestring}': "
            "only Life-like rules without B0 or B1 are supported."
        )
    birth_digits, survival_digits = match.groups()

    if not survival_digits:
        raise ValueError(
            f"Invalid rule '{rulestring}': "
            "no survival conditions."
        )
    def increasing(s): return all(s[i] < s[i+1] for i in range(len(s)-1))
    if not increasing(birth_digits) or not increasing(survival_digits):
        raise ValueError(
            f"Invalid rule '{rulestring}': "
            "only Life-like rules are supported."
        )

    birth_at = [False] * 9
    survival_at = [False] * 9
    for digit in birth_digits:
        birth_at[int(digit)] = True
    for digit in survival_digits:
        survival_at[int(digit)] = True

    return birth_at, survival_at


def build_ruletable(rulestring):
    birth_at, survival_at = get_rule(rulestring)
    ruletable = ([  f"@RULE Rotor{rulestring}",
                    "@TABLE",
                    "n_states:19",
                    "neighborhood:Moore",
                    "symmetries:permute"
    ])
    for i in range(8):
        ruletable.append(f"var on{i} = {{1,3,5,7,9,11,13,15,17}}")
        ruletable.append(f"var off{i} = {{0,2,4,6,8,10,12,14,16,18}}")

    def stator_neighbors(state):
        return (state - 1) // 2

    for total_neighbors in range(9):
        for state in range(1,19):
            rotor_neighbors = total_neighbors - stator_neighbors(state)
            if rotor_neighbors < 0:
                continue

            cell_is_on = (state % 2 == 1)
            if ( (cell_is_on and not survival_at[total_neighbors])
                 or (not cell_is_on and birth_at[total_neighbors])
            ):
                target = state + 1 if cell_is_on else state - 1
                neighbors = ",".join([
                    f"on{i}" if i < rotor_neighbors
                    else f"off{i}" for i in range(8)
                ])
                ruletable.append(f"{state}, {neighbors}, {target}")

    return '\n'.join(ruletable)
