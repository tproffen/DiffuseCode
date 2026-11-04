#!/bin/bash
date > parallel.date
mpiexec -n 5 discus_suite -macro suite.mac > LOG.LOG
date >> parallel.date
