from .rules import get_rulestring_from_pattern


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
        self.stator_boundary = lt.pattern("", "b12345678s012345678")
        self.adjacent_rotor = None      # Rotor cells with a stator neighbor.
        self.rotor_phases = None
        self.change_envelopes = None

    @property
    def initial_population(self):
        if self.initial.getrule().endswith("History"):
            return self.initial_two_state.population
        else:
            return float("inf")

    # Restrict the stator to the search area
    def make_stator(self, adjustments, distance):
        adjust_left, adjust_right, adjust_top, adjust_bottom = adjustments

        self.stator[    0 - adjust_left : self.width + adjust_right,
                        0 - adjust_top : self.height + adjust_bottom
                   ] = 1
        self.stator -= self.rotor
        self.stator_boundary += self.stator

        if distance is not None and self.rotor.nonempty():
            distance_mask = (
                self.stator.owner.pattern("", "b1e2-cn3-c4-c5678s012345678")
            )
            distance_mask += self.rotor
            self.stator &= distance_mask[distance]
            if distance == 0:
                self.stator_boundary = self.stator_boundary[1] & self.rotor
            else:
                self.stator_boundary &= distance_mask[distance]

        self.stator_boundary = (self.stator_boundary[1]
                                - self.stator_boundary
                                - self.rotor)
        self.stator += self.stator_boundary

        self.adjacent_rotor = self.stator[1] & (self.rotor - self.rotor[1])

    def analyze_pattern(self, ticks, adjustments, distance):
        pattern_is_history = self.initial.getrule().endswith("History")
        if pattern_is_history:
            final_pattern = self.initial[ticks]
            envelope = flatten_multistate(final_pattern)
            self.rotor += flatten_multistate(final_pattern ^ self.initial)
            self.initial_stator_on += envelope ^ self.rotor
        else:
            self.rotor += flatten_multistate(self.initial)

        self.make_stator(adjustments, distance)

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
