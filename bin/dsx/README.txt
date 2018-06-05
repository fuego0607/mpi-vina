
#============================================================================
# DSX; Knowledge-based scoring function for the assessment
#      of receptor-ligand interactions
# Copyright (C) 2009, 2010, 2011  Gerd Neudert and Gerhard Klebe
#
# HotspotsX; Generating acnt contour maps based on DSX pair potentials
# Copyright (C) 2009, 2010, 2011  Gerd Neudert and Gerhard Klebe
#
# Usage of DSX and HotspotsX is free without any limitations.
# Redistribution is NOT allowed without explicit permission.
# The programs are distributed in the hope that they will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#----------------------------------------------------------------------------
#
# author:      Gerd Neudert
#              gneudert(place_at_here)web.de
# supervisor:  Prof. Dr. Gerhard Klebe
#              klebe(place_at_here)staff.uni-marburg.de
# Department of Pharmaceutical Chemistry, Philipps-University of Marburg,
# Marbacher Weg 6, 35037 Marburg, Germany
#---------------------------------------------------------------------------- 


Installation notes:
===================

 -> If you read this, you have successfully downloaded, and uncompressed
     the necessary files:
        -> Precompiled DSX executables for Linux 32Bit, Linux 64Bit
            OR Intel-based Mac 32Bit, 64BIT for both, DSX and HotspotsX.
        -> A Directory containing PDB pair- and SR-potentials.
        -> An example maps file for hotspotsX: ACC_DON_AnD_HYD_ARO_map.def

 -> Please note, that we have tested the programs only on linux platforms
     (Field reports from Mac users are welcome!)

 -> Execute the program with option '-h' to get detailed and hopefully
     self-explaining usage information.

 -> Detailed explanation of the applied potentials and their theoretical
     background will be published soon.

 -> If you like to use the CSD-based pair-, torsion- and SR-potentials,
     there are two options:
        -> You can use them via our web-interface DSX-online:
            http://www.agklebe.de/drugscore
        -> After publication of DSX you can download them from the
            Cambridge Crystallographic Data Centre (CCDC), if you
            are holding a valid CSD license. We will provide a link
            as soon as the publication is out.

 -> Dependencies:
        -> Linux: 
           NOTE: If you have problems with library versions, you will
                 find compatible versions in the same directory where
                 the binary resides (set your LD_LIBRARY_PATH to this
                 directory).
                 Due to some symbol looup problems I used gcc4.3 for
                 32bit versions. However, for usage of 64bit versions
                 you need gcc4.4 or higher.
                 If even the 32bit version does not work for you,
                 please try the RHEL version.
            - libstdc++.so.6
            - libm.so.6
            - libgcc_s.so.1
            - libc.so.6
            - libpthread.so.0 (only for HotspotsX)
        -> RHEL Linux 5 and higher:
           As some people still have problems running the standard
           Linux version, I added an additional version compiled
           with gcc4.1.2 under Scientific Linux 5.6.
        -> MacOS:
            - /usr/lib/libstdc++.6.dylib (compatibility version 7.0.0)
            - /usr/lib/libSystem.B.dylib (compatibility version 1.0.0)


Marburg 16.08.2011,
Gerd Neudert

