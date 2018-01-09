// -----------------------------------------------------------------------------
//  Copyright (c) 2013-2018 Sun Yat-Sen University (SYSU). All Rights Reserved.
//
//  SYSU grants permission to use, copy, modify, and distribute this software
//  and its documentation for NON-COMMERCIAL purposes and without fee, provided
//  that this copyright notice appears in all copies.
//
//  SYSU provides this software "as is," without representations or warranties
//  of any kind, either expressed or implied, including but not limited to the
//  implied warranties of merchantability, fitness for a particular purpose, 
//  and noninfringement. SYSU shall not be liable for any damages arising from
//  any use of this software.
//
//  Authors: Qiang Huang (huangq2011@gmail.com)
//           Guihong Ma  (maguihong@vip.qq.com)
//
//  Created on:       2016-01-02
//  Last Modified on: 2018-01-06
//  Version 1.0.0
// -----------------------------------------------------------------------------

#include <stdio.h>
#include <errno.h>

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cfloat>
#include <cmath>
#include <ctime>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

// -----------------------------------------------------------------------------
//  For Linux directory
// -----------------------------------------------------------------------------
#include <unistd.h>
#include <stdarg.h>
#include <sys/stat.h>
#include <sys/types.h>
#define LINUX_

// -----------------------------------------------------------------------------
//  For Windows directory
// -----------------------------------------------------------------------------
// #include <direct.h>
// #include <io.h>

#include "def.h"
#include "util.h"
#include "pri_queue.h"
#include "qalsh_col.h"
#include "l2_alsh.h"
#include "l2_alsh2.h"
#include "xbox.h"
#include "h2_alsh.h"
#include "srp_lsh.h"
#include "sign_alsh.h"
#include "simple_lsh.h"
#include "amip.h"

using namespace std;
