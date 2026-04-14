//
// Created by themi on 10/4/2026.
//

/*
 * Implemenataion based on:
 * A. Danalis, G. Marin, C. McCurdy, J. S. Meredith, P. C. Roth, K. Spafford,
 * V. Tipparaju, and J. S. Vetter.
 * The scalable heterogeneous computing (shoc) benchmark suite.
 * In Proceedings of the 3rd Workshop on General-Purpose Computation on
 * Graphics Processing Units, 2010.
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include "md_pipelined.h"

int populateNeighborList(TYPE currDist[maxNeighbors],
                         int currList[maxNeighbors],
                         const int i,
                         int NL[nAtoms][maxNeighbors]);

int buildNeighborList(TYPE position_x[nAtoms],
                      TYPE position_y[nAtoms],
                      TYPE position_z[nAtoms],
                      int NL[nAtoms][maxNeighbors]);

/* -----------------------------------------------------------------------
 * distance — υπολογισμός τετραγώνου απόστασης μεταξύ ατόμων i και j
 * ----------------------------------------------------------------------- */
TYPE distance(TYPE position_x[nAtoms],
              TYPE position_y[nAtoms],
              TYPE position_z[nAtoms],
              int i, int j)
{
    TYPE delx, dely, delz;
    delx = position_x[i] - position_x[j];
    dely = position_y[i] - position_y[j];
    delz = position_z[i] - position_z[j];
    return delx*delx + dely*dely + delz*delz;
}

/* -----------------------------------------------------------------------
 * populateNeighborList — αντιγράφει την τοπική λίστα γειτόνων στον NL
 * ----------------------------------------------------------------------- */
int populateNeighborList(TYPE currDist[maxNeighbors],
                         int currList[maxNeighbors],
                         const int i,
                         int NL[nAtoms][maxNeighbors])
{
    int neighborIter, validPairs;
    validPairs = 0;
    for (neighborIter = 0; neighborIter < maxNeighbors; neighborIter++) {
        NL[i][neighborIter] = currList[neighborIter];
        validPairs++;
    }
    return validPairs;
}

/* -----------------------------------------------------------------------
 * buildNeighborList — για κάθε άτομο i βρίσκει τους maxNeighbors γείτονες
 * ----------------------------------------------------------------------- */
int buildNeighborList(TYPE position_x[nAtoms],
                      TYPE position_y[nAtoms],
                      TYPE position_z[nAtoms],
                      int NL[nAtoms][maxNeighbors])
{
    int totalPairs, i, j, k;
    TYPE distIJ;

    totalPairs = 0;

    for (i = 0; i < nAtoms; i++) {
        int  currList[maxNeighbors];
        TYPE currDist[maxNeighbors];

        for (k = 0; k < maxNeighbors; k++) {
            currList[k] = 0;
            currDist[k] = 999999999.0;
        }

        for (j = 0; j < maxNeighbors; j++) {
            if (i == j) continue;
            distIJ       = distance(position_x, position_y, position_z, i, j);
            currList[j]  = j;
            currDist[j]  = distIJ;
        }

        totalPairs += populateNeighborList(currDist, currList, i, NL);
    }
    return totalPairs;
}

/* -----------------------------------------------------------------------
 * main
 * ----------------------------------------------------------------------- */
int main(void)
{
    int i, j, iter, totalPairs;

    srand(8650341L);

    TYPE position_x[nAtoms];
    TYPE position_y[nAtoms];
    TYPE position_z[nAtoms];
    TYPE force_x[nAtoms];
    TYPE force_y[nAtoms];
    TYPE force_z[nAtoms];

    int NL[nAtoms][maxNeighbors];

    int32_t neighborList[nAtoms * maxNeighbors];

    for (i = 0; i < nAtoms; i++) {
        position_x[i] = (TYPE)rand() / RAND_MAX;
        position_y[i] = (TYPE)rand() / RAND_MAX;
        position_z[i] = (TYPE)rand() / RAND_MAX;
        force_x[i]    = 0.0;
        force_y[i]    = 0.0;
        force_z[i]    = 0.0;
    }

    for (i = 0; i < nAtoms; i++)
        for (j = 0; j < maxNeighbors; j++)
            NL[i][j] = 0;

    totalPairs = buildNeighborList(position_x, position_y, position_z, NL);
    printf("Total neighbor pairs: %d\n", totalPairs);

    for (i = 0; i < nAtoms; i++)
        for (j = 0; j < maxNeighbors; j++)
            neighborList[i * maxNeighbors + j] = (int32_t)NL[i][j];

    for (iter = 0; iter < MAX_ITERATION; iter++) {
        md_kernel(force_x, force_y, force_z,
                  position_x, position_y, position_z,
                  neighborList);
    }

    printf("=== Force results (first 4 atoms) ===\n");
    for (i = 0; i < 4; i++) {
        printf("force[%d] = (%f, %f, %f)\n",
               i, force_x[i], force_y[i], force_z[i]);
    }

    return 0;
}