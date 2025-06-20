C
C    Copyright (C) 2019 Weierstrass-Institut fuer
C                       Angewandte Analysis und Stochastik (WIAS)
C
C    Author:  Joerg Polzehl
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.
C
C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.
C
C  You should have received a copy of the GNU General Public License
C  along with this program; if not, write to the Free Software
C  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
C  USA.
C
C  The following routines were part of the aws package and contain
C  FORTRAN 77 code needed in R functions aws, vaws,
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


#include <stdio.h>
#include <math.h>
#include <stdbool.h>

// External functions
extern double lkern(int, double);
extern double KLdistsi(double*, double*, double*, int);

void pvaws(
    double *y,
    double *pos,
    int nv,
    int nvd,
    int n1,
    int n2,
    int n3,
    double hakt,
    double lambda,
    double *theta,
    double *bi,
    double *bin,
    double *thnew,
    double *invcov,
    int ncores,
    double spmin,
    double *lwght,
    double *wght,
    double *swjy,
    int np1,
    int np2,
    int np3
) {
    int ih1, ih2, ih3, i1, i2, i3, j1, j2, j3, jw1, jw2, jw3, jwind3, jwind2;
    int iind, jind, jind3, jind2, clw1, clw2, clw3, dlw1, dlw2, dlw3, dlw12, n12, k, thrednr, iindp, jindp, ipindp, jpindp;
    double sij, swj, z1, z2, z3, wj, hakt2, w1, w2, sijp;
    int np1, np2, np3;
    int ip1, ip2, ip3, nph1, nph2, nph3, ipind, jp1, jp2, jp3, jpind;
    bool aws;
    double spf;

    thrednr = 1;
    hakt2 = hakt * hakt;
    spf = 1.0 / (1.0 - spmin);
    ih1 = floor(hakt);
    aws = lambda < 1e35;

    // First calculate location weights
    w1 = wght[0];
    w2 = wght[1];
    ih3 = floor(hakt / w2);
    ih2 = floor(hakt / w1);
    ih1 = floor(hakt);
    if (n3 == 1) ih3 = 0;
    if (n2 == 1) ih2 = 0;
    clw1 = ih1;
    clw2 = ih2;
    clw3 = ih3;
    dlw1 = ih1 + clw1 + 1;
    dlw2 = ih2 + clw2 + 1;
    dlw3 = ih3 + clw3 + 1;
    dlw12 = dlw1 * dlw2;
    nph1 = (np1 - 1) / 2;
    nph2 = (np2 - 1) / 2;
    nph3 = (np3 - 1) / 2;
    n12 = n1 * n2;
    z2 = 0.0;
    z3 = 0.0;

    for (j3 = -clw3; j3 <= clw3; j3++) {
        if (n3 > 1) {
            z3 = j3 * w2;
            z3 = z3 * z3;
            ih2 = floor(sqrt(hakt2 - z3) / w1);
            jind3 = (j3 + clw3) * dlw12;
        } else {
            jind3 = 0;
        }
        for (j2 = -ih2; j2 <= ih2; j2++) {
            if (n2 > 1) {
                z2 = j2 * w1;
                z2 = z3 + z2 * z2;
                ih1 = floor(sqrt(hakt2 - z2));
                jind2 = jind3 + (j2 + clw2) * dlw1;
            } else {
                jind2 = 0;
            }
            for (j1 = -ih1; j1 <= ih1; j1++) {
                jind = j1 + clw1 + 1 + jind2;
                z1 = j1;
                lwght[jind] = lkern(2, (z1 * z1 + z2) / hakt2);
            }
        }
    }

    // Rescale bi with 1/lambda
    for (iind = 0; iind < n1 * n2 * n3; iind++) {
        iindp = pos[iind];
        if (iindp == 0) continue;
        bi[iindp] = bi[iindp] / lambda;
    }

    // Call to rchkusr() - assuming it's a placeholder for some user-defined function
    // rchkusr();

    for (iind = 0; iind < n1 * n2 * n3; iind++) {
        iindp = pos[iind];
        if (iindp == 0) continue;

        i1 = iind % n1;
        if (i1 == 0) i1 = n1;
        i2 = (iind - i1) / n1 + 1;
        if (i2 == 0) i2 = n2;
        i3 = (iind - i1 - (i2 - 1) * n1) / n12 + 1;

        swj = 0.0;
        for (k = 0; k < nv; k++) {
            swjy[k * ncores + thrednr] = 0.0;
        }

        for (jw3 = -clw3; jw3 <= clw3; jw3++) {
            j3 = jw3 + i3;
            if (j3 < 1 || j3 > n3) continue;
            jwind3 = (jw3 + clw3) * dlw12;
            jind3 = (j3 - 1) * n12;
            z3 = jw3 * w2;
            z3 = z3 * z3;
            if (n2 > 1) ih2 = floor(sqrt(hakt2 - z3) / w1);
            for (jw2 = -ih2; jw2 <= ih2; jw2++) {
                j2 = jw2 + i2;
                if (j2 < 1 || j2 > n2) continue;
                jwind2 = jwind3 + (jw2 + clw2) * dlw1;
                jind2 = (j2 - 1) * n1 + jind3;
                z2 = jw2 * w1;
                z2 = z3 + z2 * z2;
                ih1 = floor(sqrt(hakt2 - z2));
                for (jw1 = -ih1; jw1 <= ih1; jw1++) {
                    j1 = jw1 + i1;
                    if (j1 < 1 || j1 > n1) continue;
                    jind = j1 + jind2;
                    jindp = pos[jind];
                    if (jindp == 0) continue;
                    wj = lwght[jw1 + clw1 + 1 + jwind2];
                    if (aws) {
                        sij = 0.0;
                        for (ip1 = i1 - nph1; ip1 <= i1 + nph1; ip1++) {
                            if (ip1 <= 0 || ip1 > n1) continue;
                            jp1 = ip1 + jw1;
                            if (jp1 <= 0 || jp1 > n1) continue;
                            for (ip2 = i2 - nph2; ip2 <= i2 + nph2; ip2++) {
                                if (ip2 <= 0 || ip2 > n2) continue;
                                jp2 = ip2 + jw2;
                                if (jp2 <= 0 || jp2 > n2) continue;
                                for (ip3 = i3 - nph3; ip3 <= i3 + nph3; ip3++) {
                                    if (sij > 1.0) continue;
                                    if (ip3 <= 0 || ip3 > n3) continue;
                                    ipind = ip1 + (ip2 - 1) * n1 + (ip3 - 1) * n12;
                                    ipindp = pos[ipind];
                                    if (ipindp == 0) continue;
                                    jp3 = ip3 + jw3;
                                    if (jp3 <= 0 || jp3 > n3) continue;
                                    jpind = jp1 + (jp2 - 1) * n1 + (jp3 - 1) * n12;
                                    jpindp = pos[jpind];
                                    if (jpindp == 0) continue;
                                    sijp = KLdistsi(&theta[jpindp * nv], &theta[ipindp * nv], &invcov[ipindp * nvd], nv);
                                    sij = fmax(sij, bi[ipindp] * sijp);
                                }
                            }
                        }
                        if (sij >= 1.0) continue;
                        if (sij > spmin) wj = wj * (1.0 - spf * (sij - spmin));
                    }
                    swj += wj;
                    for (k = 0; k < nv; k++) {
                        swjy[k * ncores + thrednr] += wj * y[k * ncores + jindp];
                    }
                }
            }
        }
        for (k = 0; k < nv; k++) {
            thnew[k * ncores + iindp] = swjy[k * ncores + thrednr] / swj;
        }
        bin[iindp] = swj;
    }
}


void pvawsme(double *y, double *yd, double *pos, int nv, int nvd, int nd, int n1, int n2, int n3, 
             double hakt, double lambda, double *theta, double *bi, double *bin, double *thnew, 
             double *ydnew, double *invcov, int ncores, double spmin, double *lwght, double *wght, 
             double *swjy, double *swjd, int np1, int np2, int np3) {

    int ih1, ih2, ih3, i1, i2, i3, j1, j2, j3, jw1, jw2, jw3, jwind3, jwind2;
    int iind, jind, jind3, jind2, clw1, clw2, clw3, dlw1, dlw2, dlw3, dlw12, n12, k, thrednr, iindp, jindp, ipindp, jpindp;
    double sij, swj, z1, z2, z3, wj, hakt2, w1, w2, sijp;
    int np1, np2, np3;
    int ip1, ip2, ip3, nph1, nph2, nph3, ipind, jp1, jp2, jp3, jpind;
    bool aws;
    double spf;

    thrednr = 1;
    hakt2 = hakt * hakt;
    spf = 1.0 / (1.0 - spmin);
    ih1 = floor(hakt);
    aws = lambda < 1e35;

    // First calculate location weights
    w1 = wght[0];
    w2 = wght[1];
    ih3 = floor(hakt / w2);
    ih2 = floor(hakt / w1);
    ih1 = floor(hakt);
    if (n3 == 1) ih3 = 0;
    if (n2 == 1) ih2 = 0;
    clw1 = ih1;
    clw2 = ih2;
    clw3 = ih3;
    dlw1 = ih1 + clw1 + 1;
    dlw2 = ih2 + clw2 + 1;
    dlw3 = ih3 + clw3 + 1;
    dlw12 = dlw1 * dlw2;
    nph1 = (np1 - 1) / 2;
    nph2 = (np2 - 1) / 2;
    nph3 = (np3 - 1) / 2;
    n12 = n1 * n2;
    z2 = 0.0;
    z3 = 0.0;

    for (j3 = -clw3; j3 <= clw3; j3++) {
        if (n3 > 1) {
            z3 = j3 * w2;
            z3 = z3 * z3;
            ih2 = floor(sqrt(hakt2 - z3) / w1);
            jind3 = (j3 + clw3) * dlw12;
        } else {
            jind3 = 0;
        }
        for (j2 = -ih2; j2 <= ih2; j2++) {
            if (n2 > 1) {
                z2 = j2 * w1;
                z2 = z3 + z2 * z2;
                ih1 = floor(sqrt(hakt2 - z2));
                jind2 = jind3 + (j2 + clw2) * dlw1;
            } else {
                jind2 = 0;
            }
            for (j1 = -ih1; j1 <= ih1; j1++) {
                jind = j1 + clw1 + 1 + jind2;
                z1 = j1;
                lwght[jind] = lkern(2, (z1 * z1 + z2) / hakt2);
            }
        }
    }

    // Rescale bi with 1/lambda
    for (iind = 0; iind < n1 * n2 * n3; iind++) {
        iindp = pos[iind];
        if (iindp == 0) continue;
        bi[iindp] = bi[iindp] / lambda;
    }

    // Call to rchkusr() - assuming it's a placeholder for some user-defined function
    // rchkusr();

    for (iind = 0; iind < n1 * n2 * n3; iind++) {
        iindp = pos[iind];
        if (iindp == 0) continue;

        i1 = iind % n1;
        if (i1 == 0) i1 = n1;
        i2 = (iind - i1) / n1 + 1;
        if (i2 == 0) i2 = n2;
        i3 = (iind - i1 - (i2 - 1) * n1) / n12 + 1;

        swj = 0.0;
        for (k = 0; k < nv; k++) {
            swjy[k * ncores + thrednr] = 0.0;
        }
        for (k = 0; k < nd; k++) {
            swjd[k * ncores + thrednr] = 0.0;
        }

        for (jw3 = -clw3; jw3 <= clw3; jw3++) {
            j3 = jw3 + i3;
            if (j3 < 1 || j3 > n3) continue;
            jwind3 = (jw3 + clw3) * dlw12;
            jind3 = (j3 - 1) * n12;
            z3 = jw3 * w2;
            z3 = z3 * z3;
            if (n2 > 1) ih2 = floor(sqrt(hakt2 - z3) / w1);
            for (jw2 = -ih2; jw2 <= ih2; jw2++) {
                j2 = jw2 + i2;
                if (j2 < 1 || j2 > n2) continue;
                jwind2 = jwind3 + (jw2 + clw2) * dlw1;
                jind2 = (j2 - 1) * n1 + jind3;
                z2 = jw2 * w1;
                z2 = z3 + z2 * z2;
                ih1 = floor(sqrt(hakt2 - z2));
                for (jw1 = -ih1; jw1 <= ih1; jw1++) {
                    j1 = jw1 + i1;
                    if (j1 < 1 || j1 > n1) continue;
                    jind = j1 + jind2;
                    jindp = pos[jind];
                    if (jindp == 0) continue;
                    wj = lwght[jw1 + clw1 + 1 + jwind2];
                    if (aws) {
                        sij = 0.0;
                        for (ip1 = i1 - nph1; ip1 <= i1 + nph1; ip1++) {
                            if (ip1 <= 0 || ip1 > n1) continue;
                            jp1 = ip1 + jw1;
                            if (jp1 <= 0 || jp1 > n1) continue;
                            for (ip2 = i2 - nph2; ip2 <= i2 + nph2; ip2++) {
                                if (ip2 <= 0 || ip2 > n2) continue;
                                jp2 = ip2 + jw2;
                                if (jp2 <= 0 || jp2 > n2) continue;
                                for (ip3 = i3 - nph3; ip3 <= i3 + nph3; ip3++) {
                                    if (sij > 1.0) continue;
                                    if (ip3 <= 0 || ip3 > n3) continue;
                                    ipind = ip1 + (ip2 - 1) * n1 + (ip3 - 1) * n12;
                                    ipindp = pos[ipind];
                                    if (ipindp == 0) continue;
                                    jp3 = ip3 + jw3;
                                    if (jp3 <= 0 || jp3 > n3) continue;
                                    jpind = jp1 + (jp2 - 1) * n1 + (jp3 - 1) * n12;
                                    jpindp = pos[jpind];
                                    if (jpindp == 0) continue;
                                    sijp = KLdistsi(&theta[jpindp * nv], &theta[ipindp * nv], &invcov[ipindp * nvd], nv);
                                    sij = fmax(sij, bi[ipindp] * sijp);
                                }
                            }
                        }
                        if (sij >= 1.0) continue;
                        if (sij > spmin) wj = wj * (1.0 - spf * (sij - spmin));
                    }
                    swj += wj;
                    for (k = 0; k < nv; k++) {
                        swjy[k * ncores + thrednr] += wj * y[k * ncores + jindp];
                    }
                    for (k = 0; k < nd; k++) {
                        swjd[k * ncores + thrednr] += wj * yd[k * ncores + jindp];
                    }
                }
            }
        }
        for (k = 0; k < nv; k++) {
            thnew[k * ncores + iindp] = swjy[k * ncores + thrednr] / swj;
        }
        for (k = 0; k < nd; k++) {
            ydnew[k * ncores + iindp] = swjd[k * ncores + thrednr] / swj;
        }
        bin[iindp] = swj;
    }
}


double KLdistsr(double *thi, double *thj, double *si2, int nv) {
    double z = 0.0;
    double zdk;
    int k, l, m = 0;

    for (k = 0; k < nv; k++) {
        zdk = thi[k] - thj[k];
        if (k > 0) {
            for (l = 0; l < k; l++) {
                z += 2.0 * (thi[l] - thj[l]) * zdk * si2[m];
                m++;
            }
        }
        z += zdk * zdk * si2[m];
        m++;
    }

    return z;
}

double KLdistsi(double *thi, double *thj, double *si2, int nv) {
    double z = 0.0;
    double zdk;
    int k, l, m = 0;

    for (k = 0; k < nv; k++) {
        zdk = thi[k] - thj[k];
        if (k > 0) {
            for (l = 0; l < k; l++) {
                z += 2.0 * (thi[l] - thj[l]) * zdk * si2[m];
                m++;
            }
        }
        z += zdk * zdk * si2[m];
        m++;
    }

    return z;
}

