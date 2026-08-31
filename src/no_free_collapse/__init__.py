"""No Free Collapse: exact and approximate interaction-order utilities."""

from .born import born_score, normalized_born_probability, polynomial_state, preparation_order
from .completion import (
    completion_dual_value,
    completion_matrix,
    completion_primal_trace,
    offdiag_cube_extrema,
    offdiag_cube_max,
    six_variable_range_hafnian_bound,
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
from .projection_gradient import (
    SIX_VARIABLE_EDGES,
    complementary_four_hafnians,
    rank_three_projection_gradient_defect,
    rank_three_projection_involution,
    rank_three_projection_stability_terms,
    six_variable_edge_vector,
    six_variable_gradient_energies,
    six_variable_perfect_matching_operator,
)
from .rank_two_dual import (
    evaluate_degree_four_bernstein,
    rank_two_high_pair_bernstein_coefficients,
    rank_two_quadratic_dual_certificate,
)
from .rank_two_projection import (
    rank_two_harmonic_gradient_defect,
    rank_two_harmonic_invariants,
    rank_two_heavy_bernstein_coefficients,
    rank_two_heavy_coordinate_decomposition,
    rank_two_heavy_coordinate_projection,
    rank_two_plucker_weights,
    rank_two_projection_direct_defect,
    rank_two_projection_gradient_defect,
    rank_two_stability_decomposition,
    rank_two_tensor_q2_formula,
)

__all__ = [
    "born_score", "normalized_born_probability", "polynomial_state", "preparation_order",
    "completion_dual_value", "completion_matrix", "completion_primal_trace",
    "offdiag_cube_extrema", "offdiag_cube_max", "six_variable_range_hafnian_bound",
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
    "SIX_VARIABLE_EDGES", "complementary_four_hafnians",
    "rank_three_projection_gradient_defect", "rank_three_projection_involution",
    "rank_three_projection_stability_terms", "six_variable_edge_vector",
    "six_variable_gradient_energies", "six_variable_perfect_matching_operator",
    "evaluate_degree_four_bernstein", "rank_two_high_pair_bernstein_coefficients",
    "rank_two_quadratic_dual_certificate",
    "rank_two_harmonic_gradient_defect", "rank_two_harmonic_invariants",
    "rank_two_heavy_bernstein_coefficients", "rank_two_heavy_coordinate_decomposition",
    "rank_two_heavy_coordinate_projection", "rank_two_plucker_weights",
    "rank_two_projection_direct_defect", "rank_two_projection_gradient_defect",
    "rank_two_stability_decomposition", "rank_two_tensor_q2_formula",
]
