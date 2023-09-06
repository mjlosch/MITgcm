## Various tests of the Weddell Sea configuration for SOCHIC

- run00: tests

- run04: "reference" forward run for 62 years (1958/01/01-2019/12/31),
  JRA55, nonlinear free surface, r-star, etc. not possible to use all
  of that in AD runs

- adrun00: short test with 20cpu
- adrun01: 10year integration with mxlMaxFlag = 3 (works also in AD)
           with full JRA55 data (too slow)
- adrun02: 10year integration with mxlMaxFlag = 1 and coarse JRA55 data
- adrun03: like adrun02 with mxlMaxFlag = 2
- adrun04: like adrun02 with mxlMaxFlag = 3
- adrun05: performance test for disk access, uses mxlMaxFlag = 1
- adrun05a: same as adrun05, use 40 instead of 20 cpu