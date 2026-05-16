# snow_like_affine_approximation

This repository contains the auxiliary search code and raw supporting files for the paper **"Best Affine Approximations on SNOW-Like Ciphers"**.

The files are intended to support the affine-approximation search and the signed bias aggregation reported for SNOW 2.0 and SNOW 3G. The paper reports all final comparison values in the **bias** convention. If a cited work reports a correlation `c`, the corresponding bias is `|c|/2`.

## Requirements

The search programs are C++ programs using the Gurobi C++ API. Re-running the searches requires:

- a C++ compiler compatible with the local Gurobi installation;
- a valid Gurobi license;
- the Gurobi C++ headers and libraries configured in the local build environment.

No fixed build command is provided here because Gurobi include and library paths are installation-dependent. When re-running a program, use the corresponding cipher directory as the working directory so that the program can read the inequality and output files in that directory.

## Directory Structure

```text
snow_like_affine_approximation-main/
├── README.md
├── snow 2.0/
└── snow 3g/
```

## SNOW 2.0 Files

Directory: `snow 2.0/`

- `main.cpp`: C++ implementation of the MILP-based search used for the SNOW 2.0 FSM approximation.
- `equation.h`: helper definitions used by `main.cpp`.
- `ineq_aes_25.txt`: linear-inequality description for the AES S-box LAT model used in the MILP formulation.
- `all_objvalue1.txt`: raw objective-value output from the search stage.
- `all_sol1.txt`: raw solution-pool mask assignments corresponding to the search output.
- `sol_gather1.txt`: the 40 non-zero contributing mask entries used for the signed aggregation under the selected fixed masks.
- `聚集结果.xlsx`: spreadsheet record used to check the SNOW 2.0 aggregation bookkeeping.

In the paper, the SNOW 2.0 result is the aggregated bias `2^{-15.03}` for the selected fixed masks. The objective values in `all_objvalue1.txt` are internal search scores and should not be read directly as the final aggregated bias. The final bias is obtained only after summing the signed trail contributions for the selected fixed masks.

## SNOW 3G Files

Directory: `snow 3g/`

- `main.cpp`: C++ implementation of the MILP-based search used for the SNOW 3G FSM approximation.
- `equation.h`: helper definitions used by `main.cpp`.
- `ineq_aes_25.txt`: linear-inequality description for the AES S-box LAT model.
- `ineq_sq_21.txt`: linear-inequality description for the SNOW 3G `S_Q` S-box LAT model.
- `ineq_xor_11.txt` and `ineq_xor_5.txt`: linear-inequality descriptions for XOR constraints used by the MILP model.
- `mask.txt`: mask data used by the SNOW 3G search/aggregation workflow, stored in the program's row-wise binary format.
- `obj_gather.txt`: objective values for the 100 paths collected for aggregation.
- `all_sol_gather.txt`: mask assignments for the same 100 paths, stored in the program's row-wise binary format.

In the paper, the SNOW 3G result is the aggregated bias `2^{-21.66}` for the binary FSM approximation over three consecutive keystream words. The appendix lists only the 20 largest contributing paths to keep the paper concise. The complete 100-path support is provided by `obj_gather.txt` and `all_sol_gather.txt`.

As in the SNOW 2.0 case, the objective values in `obj_gather.txt` are internal search scores. They support the path selection and aggregation check, but the final paper value is the signed aggregated bias, not a direct transcription of the largest objective value.

## Relation to the Paper

- SNOW 2.0: the support files correspond to the search and aggregation used for the bias `2^{-15.03}`.
- SNOW 3G: the support files correspond to the search and aggregation used for the bias `2^{-21.66}`.
- The fast-correlation attack complexities reported in the paper are computed separately from the bias search, using published 4-tree estimation formulas. They are not output directly by the C++ search programs in this repository.

## Notes for Verification

1. Check the fixed-mask relation and the contributing paths in the paper first.
2. Use the raw solution and objective files in the corresponding cipher directory to verify that the listed paths are included in the support data.
3. Use the signed aggregation convention described in the paper when converting contributing paths into the final bias.
4. Use the paper's fast-correlation formulas, not the MILP objective values, to verify the reported data/time/memory exponents.
