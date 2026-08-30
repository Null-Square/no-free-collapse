"""No Free Collapse: exact and approximate interaction-order utilities."""

from .born import born_score, normalized_born_probability, polynomial_state, preparation_order
from .completion import (
    completion_dual_value,
    completion_matrix,
    completion_primal_trace,
    offdiag_cube_max,
)
from .conditioning import (
    centered_relative_variation,
    condition_number,
    delta_from_condition_number,
    required_condition_number_for_coefficient,
    spectral_coefficient_bound,
    spectral_tail_energy_bound,
    truncation_bound,
)
from .constructions import parity_state_coefficients
from .gram import (
    effect_gram,
    feature_vector,
    gram_probability,
    monomial_masks,
    optimal_absolute_linear_readout_for_fixed_Q,
    optimal_linear_readout_for_fixed_Q,
    realize_gram_pair,
)
from .hafnian_bounds import (
    bounded_quadratic_hafnian_bound,
    disjoint_pair_basis_blocks,
    disjoint_pair_hafnian_directional_derivative,
    disjoint_pair_hafnian_value,
    disjoint_pair_quadratic,
    disjoint_pair_tangent_invariants,
    full_parity_power_coefficient,
    hafnian,
    six_variable_critical_hafnian_quadratic_coefficient,
    six_variable_critical_path_second_order_bound,
    twice_offdiag,
)
from .interactions import interaction_degree, walsh_coefficient, walsh_coefficients
from .paired import matched_pair_asymptotic_rate, matched_pair_capacity, matched_pair_gram

__all__ = [
    "born_score", "normalized_born_probability", "polynomial_state", "preparation_order",
    "completion_dual_value", "completion_matrix", "completion_primal_trace", "offdiag_cube_max",
    "centered_relative_variation", "condition_number", "delta_from_condition_number",
    "required_condition_number_for_coefficient", "spectral_coefficient_bound",
    "spectral_tail_energy_bound", "truncation_bound", "parity_state_coefficients",
    "effect_gram", "feature_vector", "gram_probability", "monomial_masks",
    "optimal_absolute_linear_readout_for_fixed_Q", "optimal_linear_readout_for_fixed_Q",
    "realize_gram_pair", "bounded_quadratic_hafnian_bound", "disjoint_pair_basis_blocks",
    "disjoint_pair_hafnian_directional_derivative", "disjoint_pair_hafnian_value",
    "disjoint_pair_quadratic", "disjoint_pair_tangent_invariants",
    "full_parity_power_coefficient", "hafnian",
    "six_variable_critical_hafnian_quadratic_coefficient",
    "six_variable_critical_path_second_order_bound", "twice_offdiag",
    "interaction_degree", "walsh_coefficient", "walsh_coefficients",
    "matched_pair_asymptotic_rate", "matched_pair_capacity", "matched_pair_gram",
]
