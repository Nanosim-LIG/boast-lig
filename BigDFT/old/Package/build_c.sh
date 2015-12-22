#!/bin/sh -x
icc -I/usr/include -ansi-alias -std=c99 -openmp -vec-report=3  -O3 -ipo -DNDEBUG -DAVX3 -mP2OPT_hlo_prefetch=F -mmic -mGLOB_default_function_attrs=\"use_fast_math=on\" -o benchmark_c.exe benchmark_c.c convolut_kinetic_per_c_k.c barrier.c
