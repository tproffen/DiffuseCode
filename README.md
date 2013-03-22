Diffuse program collection
----------------------------------------------------------------

Welcome to the 'diffuse program' collection. Depending on the
package you have downloaded you will find one or all of the
following program directories after unpacking the archive.

DISCUS  : Diffuse Scattering & Defect Structure Simulation
AUTHORS : R.B. Neder  (reinhard.neder@krist.uni-erlangen.de)
          Th. Proffen (tproffen@ornl.gov)

DIFFEV  : Generic refinement program using evolutionary algoritm
AUTHOR  : R.B. Neder  (reinhard.neder@krist.uni-erlangen.de)

KUPLOT  : General plotting program (well suited for DISCUS outp)
AUTHOR  : Th. Proffen (tproffen@ornl.gov)
          R.B. Neder  (reinhard.neder@krist.uni-erlangen.de)

MIXSCAT : Program to generate differential PDFs from n/X data
AUTHOR  : C. Wurden K. Page A. Llobet
          Th. Proffen (tproffen@ornl.gov)

INSTALLATION
----------------------------------------------------------------

To build from the source code, you need gfortran (> 4.7.x) as
well as cmake installed. Here is the simple set of commands 
to build the programs from the source:

Download the source code from discus.sourceforge.net or fork
our repository from GitHub https://github.com/tproffen/DiffuseCode

 Goto working directory and unpack

>  cd your-working-directory
>  tar -xvzf Diffuse-source-xxxxxx.tar.gz

Create build directory

>  mkdir DiffuseBuild
>  cd DiffuseBuild

Invloke cmake and change entries as needed. Then press 'c'
to configure and 'g' to create the Makefiles. If you use the GUI
version of cmake use the correpoding buttons.
  
>  cmake ../DiffuseCode

Now build and install

>  make
>  sudo make install


Binary distributions as well as the source code can be found 
at http://discus.sourceforge.net.

If you like to participate in code development, join us at


