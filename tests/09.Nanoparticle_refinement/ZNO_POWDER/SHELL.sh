#!/bin/bash
date > parallel.date
mpiexec -np 4 discus_suite -macro refine.mac > LOG.LOG
date >> parallel.date
