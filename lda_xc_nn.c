/*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_LDA_XC_NN   220   /* LDA NN XC Functional */

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
  NULL, NULL, 
  work_lda_nn, NULL, NULL
};
