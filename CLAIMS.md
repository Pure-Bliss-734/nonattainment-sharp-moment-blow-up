# Finite checks and paper locations

- I §2, App A: word relations, compactness, positive blocks and primal/dual witnesses. Output directory `12/`.
- I §§2–3: finite completion and singular nonclosure. Output directory `13/`.
- I §§2–3: exposed face and separation. Output directory `14/`.
- I §3, App B: isolating boxes, exact extrema and nonattainment. Output directory `15/`.
- I §3: independent root, rank and curvature checks. Output directory `16/`.
- I §§3,5: full-fiber neighborhood estimates. Output directory `17/`.
- I §§4–5: moment lift and sharp asymptotic order. Output directory `18/`.
- I §5, App C: optimal prefactors and continuous feasible families. Output directory `19/`.
- I §5: independent prefactor and interval checks. Output directory `20/`.
- I §4: exact examples and realizable Schur variations. Output directory `21/`.
- I §6: quadratic paths and the upper remainder. Output directory `22/`.
- I §6: universal lower remainder. Output directory `23/`.
- I §6, App D: constant-term strip and polynomial structure. Output directory `24/`.

Appendix B fixes the SHA-256 values of the two KKT input files and symmetry data in `stages/step15/inputs/`. These three files are distributed verbatim. The endpoint root intervals are in `15/summary.json`; the optimal prefactor enclosures are in `19/lower_prefactor.json` and `19/upper_prefactor.json`; the four constant-term bounds are in `24/summary.json`. The integral and rational finite checks support the analytical arguments on nonattainment and algebraicity.

Table 1 uses displayed decimal approximations. The upper prefactor enclosure begins 3.22362336809396584076574 times 10^-6. Its displayed mantissa 3.2236233680939658408 is rounded in the last place. The checker compares both exact endpoints at the displayed precision, and compares the four rational constant-term bounds exactly. The full enclosures are present in the reconstructed records.

## Appendix material index

The following objects are named in the paper. Paths beginning with a stage number are relative to the output of `--results`. Regenerated mathematical records are checked against frozen references before any downstream stage uses them.

- Distributed input or implementation: `stages/step15/inputs/symmetry.json`.
- Distributed input or implementation: `stages/step15/inputs/kkt_lower.json`.
- Distributed input or implementation: `stages/step15/inputs/kkt_upper.json`.
- Distributed input or implementation: `stages/step15/inputs/face_rows_lower.json`.
- Distributed input or implementation: `stages/step15/inputs/face_rows_upper.json`.
- Distributed input or implementation: `stages/step12/inputs/primal_lower.json`.
- Distributed input or implementation: `stages/step12/inputs/primal_upper.json`.
- Reconstructed record: `15/lower_root.json`.
- Reconstructed record: `15/upper_root.json`.
- Reconstructed record: `15/lower_endpoint.json`.
- Reconstructed record: `15/upper_endpoint.json`.
- Distributed input or implementation: `stages/step15/verify_step15.py`.
- Distributed input or implementation: `stages/step19/prefactor_model.py`.
- Distributed input or implementation: `stages/step19/inputs/parameters.json`.
- Reconstructed record: `19/lower_prefactor.json`.
- Reconstructed record: `19/upper_prefactor.json`.
- Reconstructed record: `19/lower_continuous_family.json`.
- Reconstructed record: `19/upper_continuous_family.json`.
- Reconstructed record: `19/formal_identities.json`.
- Reconstructed record: `19/lower_full_samples.json`.
- Reconstructed record: `19/upper_full_samples.json`.
- Distributed input or implementation: `stages/step19/verify_step19.py`.
- Distributed input or implementation: `stages/step24/inputs/parameters.json`.
- Distributed input or implementation: `stages/step24/partner_model.py`.
- Reconstructed record: `24/lower_partner_coefficients.json`.
- Reconstructed record: `24/upper_partner_coefficients.json`.
- Reconstructed record: `24/lower_finite_lower.json`.
- Reconstructed record: `24/upper_finite_lower.json`.
- Reconstructed record: `24/lower_finite_upper.json`.
- Reconstructed record: `24/upper_finite_upper.json`.
- Distributed input or implementation: `stages/step24/semialgebraic_model.py`.
- Reconstructed record: `24/model_structure.json`.
- Reconstructed record: `24/semialgebraic_formulas.json`.
- Reconstructed record: `24/examples.json`.
- Reconstructed record: `24/exact_tests.json`.
- Reconstructed record: `24/negative_controls.json`.
- Reconstructed record: `24/summary.json`.
- Distributed input or implementation: `stages/step24/verify_step24.py`.
- Distributed input or implementation: `stages/step24/stage24_sources.py`.

With a separately supplied manuscript LaTeX source, run the optional citation and cross-reference checks after reproduction:

```sh
python3 -B -S tools/check_manuscript.py manuscript.tex --output manuscript_check.json
```
