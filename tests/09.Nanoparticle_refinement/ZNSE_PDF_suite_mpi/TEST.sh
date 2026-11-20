#!/bin/bash
mpiexec --use-hwthread-cpus -np $DISCUS_NCPU discus_suite -macro refine.mac
#mpiexec --use-hwthread-cpus -np $DISCUS_NCPU discus_suite -macro refine_continue.mac
