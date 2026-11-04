#!/bin/bash
mpiexec -np 4 discus_suite -macro refine.mac
mpiexec -np 4 discus_suite -macro refine_continue.mac
