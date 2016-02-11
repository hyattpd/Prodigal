/*******************************************************************************
    PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
    Copyright (C) 2007-2016 University of Tennessee / UT-Battelle

    Code Author:  Doug Hyatt

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*******************************************************************************/

#ifndef _DPROG_H
#define _DPROG_H

#include <stdio.h>
#include <math.h>
#include "sequence.h"
#include "node.h"

#define MAX_SAM_OVLP 60
#define MAX_OPP_OVLP 200
#define MAX_NODE_DIST 500

int dprog(struct _node *, int, struct _training *, int);
void score_connection(struct _node *, int, int, struct _training *, int);
void eliminate_bad_genes(struct _node *, int, struct _training *);

#endif
