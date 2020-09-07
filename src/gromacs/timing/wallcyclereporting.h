/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#ifndef GMX_TIMING_WALLCYCLEREPORTING_H
#define GMX_TIMING_WALLCYCLEREPORTING_H

/* NOTE: None of the routines here are safe to call within an OpenMP
 * region */

#include <stdio.h>

#include <array>

#include "gromacs/utility/basedefinitions.h"

typedef struct gmx_wallcycle *gmx_wallcycle_t;
typedef struct gmx_wallclock_gpu_t gmx_wallclock_gpu_t;
struct t_commrec;

typedef std::array<double, ewcNR+ewcsNR> WallcycleCounts;
/* Convenience typedef */

WallcycleCounts wallcycle_sum(struct t_commrec *cr, gmx_wallcycle_t wc);
/* Return a vector of the sum of cycle counts over the nodes in
   cr->mpi_comm_mysim. */

void wallcycle_print(FILE *fplog, int nnodes, int npme,
                     int nth_pp, int nth_pme, double realtime,
                     gmx_wallcycle_t wc, const WallcycleCounts &cyc_sum,
                     struct gmx_wallclock_gpu_t *gpu_t);
/* Print the cycle and time accounting */

#endif
