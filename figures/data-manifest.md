# Figure data manifest

| Figure | Source data | Generator | Outputs | Status |
|---|---|---|---|---|
| HHI and recovery over four-fleet simplex | `sf_simplex_scan_results.json`, `sf_simplex_scan_summary.json` | `response_paper/figures/gen_fig_sf_simplex_scan.py` | `response_paper/figures/fig_sf_simplex_scan.pdf`, `.png`, `.svg` | Real, 200/200 certified profiles |
| Two-fleet support transition | Existing analytic path data embedded in generator | Existing generator | `response_paper/figures/fig_hhi_nonordering.pdf` | Existing |
| Mechanism-separated HHI and recovery | `sf_response_decomposition_results.json`, `sf_response_decomposition_summary.json` | `response_paper/figures/gen_fig_sf_response_decomposition.py` | `response_paper/figures/fig_sf_response_decomposition.pdf`, `.png` | Real, 600/600 certified regime-profile rows |
| Fixed-HHI spatial OD dispersion | `sf_spatial_od_results.json`, `sf_spatial_od_summary.json` | `response_paper/figures/gen_fig_sf_spatial_od.py` | `response_paper/figures/fig_sf_spatial_od.pdf`, `.png` | Real, 82/82 certified fixed-HHI profiles |
| HHI association and reversals across fleet penetration | `sf_alpha_screening_results.json`, `sf_alpha_screening_summary.json` | `response_paper/figures/gen_fig_sf_alpha_screening.py` | `response_paper/figures/fig_sf_alpha_screening.pdf`, `.png` | Real, 500/500 certified alpha-profile rows |
| Braess raw and activated response matrices | `response_paper/figures/fig_braess_response_matrices_data.json` | `response_paper/figures/gen_fig_braess_response_matrices.py` | `response_paper/figures/fig_braess_response_matrices.pdf`, `.png` | Real analytic matrix export; algebraic assertions passed |
| Sioux Falls coalition composition | `sf_coalition_composition_results.json`, `sf_coalition_composition_summary.json` | `response_paper/figures/gen_fig_sf_coalition_composition.py` | `response_paper/figures/fig_sf_coalition_composition.pdf`, `.png`, `.svg` | Real, 40/40 certified identity-preserving coalition profiles |
| Sioux Falls coalition degree | `sf_coalition_degree_results.json`, `sf_coalition_degree_summary.json` | `response_paper/figures/gen_fig_sf_coalition_degree.py` | `response_paper/figures/fig_sf_coalition_degree.pdf`, `.png`, `.svg` | Real, 54/54 certified profiles |
| Regime and cross-network validation | `regime_theorem_verification.json`, `cross_network_regime_summary.json`, `sf_coalition_design_analysis.json` | `response_paper/figures/gen_fig_regime_validation.py` | `response_paper/figures/fig_regime_validation.pdf`, `.png`, `.svg` | Real analytical and certified results; 480/480 cross-network profiles |
| Wildfire hazard regime map | Planned hazard-state summary JSON | Planned generator | best-K and regret heatmaps | Pending; no manuscript figure until certified |
| Reference-curvature transition stability | Planned V0/V1/V2 summary JSON | Planned generator | transition and cross-regret figure | Pending; no manuscript figure until certified |
| Merge-score accuracy and savings | Planned score audit JSON | Planned generator | accuracy/regret/solve-count figure | Pending; no manuscript figure until certified |

All numerical figures must be regenerated from committed machine-readable data.
No mock or planning values may enter manuscript prose.
