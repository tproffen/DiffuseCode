#!/bin/bash

echo "# Modified by Diffuse Install on `date`"         >> ${HOME}/.profile
echo "export PATH=$PATH:/usr/local/bin"                >> ${HOME}/.profile
echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib"  >> ${HOME}/.profile
echo "# End modification by Diffuse Install"           >> ${HOME}/.profile
