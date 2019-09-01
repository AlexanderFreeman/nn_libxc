/*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_LDA_XC_NN   220   /* LDA NN XC Functional */

typedef struct{
  double alpha;       /* parameters for LDA XC NN energy */
  double gamma[2];
  double beta1[2];
  double beta2[2];
  double a[2], b[2], c[2], d[2];
} lda_xc_nn_params;

static lda_xc_nn_params params_values = {
  1.0,
  {-0.1423, -0.0843},   
  { 1.0529,  1.3981}, 
  { 0.3334,  0.2611}, 
  { 0.0311,  0.01555},
  {-0.048,  -0.0269},   
  { 0.0020191519406228,  0.00069255121311694},
  {-0.0116320663789130, -0.00480126353790614}
};

static void 
lda_xc_nn_init(xc_func_type *p)
{
  lda_xc_nn_params *params;

  assert(p != NULL && p->params == NULL);
  p->params = malloc(sizeof(lda_xc_nn_params));
  params = (lda_xc_nn_params *) (p->params);

  memcpy(params, &params_values, sizeof(lda_xc_nn_params));
}

#include "maple2c/lda_xc_nn.c"

#define func maple2c_func
#include "work_lda_nn.c"


const xc_func_info_type xc_func_info_lda_xc_nn = {
  XC_LDA_XC_NN,
  XC_EXCHANGE_CORRELATION,
  "LDA XC NN Functional by A.Ryabov and P.Zhilyaev",
  XC_FAMILY_LDA,
  {&xc_ref_Rajagopal1978_L943, &xc_ref_MacDonald1979_2977, &xc_ref_Engel1995_2750, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_EXC,
  1e-24,
  0, NULL, NULL,
  lda_xc_nn_init, NULL, 
  work_lda_nn, NULL, NULL
};
