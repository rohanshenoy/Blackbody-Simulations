# Legacy Windows implementation

`bbsim1freq.py` is the original, unmodified simulation script developed on
Windows with AEDT 2023 R2 and Spyder. It is kept as the physics reference for
the `bbsim/` package: every HFSS API call in `bbsim/hfss_setup.py` and
`bbsim/extract.py` replicates a call in this file.

It is not runnable headlessly: it opens a GUI session by project name, reads
`sys.argv` at import, writes fixed scratch paths, and blocks on `input()` when
interrupted. Adaptive angular refinement (`find_regions_to_refine`,
`run_refined_plane_wave`) lives only here until it is ported and validated.

The Windows baseline settings that produced the reference `crack1Rohan` and
`crack2` datasets were `design_name = "crack2"` (or `"crack1Rohan"`),
`sweep = "discrete"`, incident theta 0 to 180 and phi 0 to 90 in 45 degree
steps, `b = 0.1` (crack2) or `0.05` (crack1Rohan), as recorded in the
`crack1Rohan_500GHz` and `crack2_500GHz` design blocks of `InfParallelPlate.aedt`.
