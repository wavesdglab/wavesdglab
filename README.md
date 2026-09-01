WavesDGlab is a set of MATLAB scripts developed to test and study
high-order continuous and discontinuous Galerkin finite elements methods.

License
-------

  WavesDGlab is copyright (C) 2023-2026, CNRS, Inria, ENSTA

  It distributed under the terms of the GNU General Public License, Version 3.

  Exception: Files in directory 'tools' are external contributions
  See LICENSE.txt for more information.

Authors
-------

  * Maintainer:
    - Axel Modave (axel.modave@ensta.fr)
  * Contributors:
    - Théophile Chaumont-Frelet (theophile.chaumont@inria.fr)
    - Pierre Marchand (pierre.marchand@inria.fr)
    - Simone Pescuma (simone.pescuma@ensta.fr)
    - Timothée Raynaud (timothee.raynaud@ensta.fr)

How to use it?
--------------

  1. Create 'setup.m' from 'setup_default.m' and define 'directoryGmsh'
  2. Run 'setup.m' in the root directory
  3. Run a script in directory 'scripts'

General structure
-----------------

  * Directory 'benchmarks': Description of the benchmarks (mesh, solution, etc.)
  * Directory 'output': Output files of the code
  * Directory 'scripts': Main routines to run
  * Directory 'sources': Finite element routines
  * Directory 'tools': External resources
