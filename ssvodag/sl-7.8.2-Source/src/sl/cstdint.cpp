//+++HDR+++
//======================================================================
//  This file is part of the SL software library.
//
//  Copyright (C) 1993-2014 by Enrico Gobbetti (gobbetti@crs4.it)
//  Copyright (C) 1996-2014 by CRS4 Visual Computing Group, Pula, Italy
//
//  For more information, visit the CRS4 Visual Computing Group 
//  web pages at http://www.crs4.it/vvr/.
//
//  This file may be used under the terms of the GNU General Public
//  License as published by the Free Software Foundation and appearing
//  in the file LICENSE included in the packaging of this file.
//
//  CRS4 reserves all rights not expressly granted herein.
//  
//  This file is provided AS IS with NO WARRANTY OF ANY KIND, 
//  INCLUDING THE WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS 
//  FOR A PARTICULAR PURPOSE.
//
//======================================================================
//---HDR---//
#include <sl/cstdint.hpp>
#include <math.h>

namespace sl {

  /**
   *  Integer Square Root function, from comp.graphics.algorithms
   *  Contributors include Arne Steinarson for the basic approximation idea, Dann 
   *  Corbit and Mathew Hendry for the first cut at the algorithm, Lawrence Kirby 
   *  for the rearrangement, improvments and range optimization and Paul Hsieh 
   *  for the round-then-adjust idea.
   */
  uint32_t isqrt(uint32_t x) {
    static const uint8_t sqq_table[] = {
      0,  16,  22,  27,  32,  35,  39,  42,  45,  48,  50,  53,  55,  57,
      59,  61,  64,  65,  67,  69,  71,  73,  75,  76,  78,  80,  81,  83,
      84,  86,  87,  89,  90,  91,  93,  94,  96,  97,  98,  99, 101, 102,
      103, 104, 106, 107, 108, 109, 110, 112, 113, 114, 115, 116, 117, 118,
      119, 120, 121, 122, 123, 124, 125, 126, 128, 128, 129, 130, 131, 132,
      133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 144, 145,
      146, 147, 148, 149, 150, 150, 151, 152, 153, 154, 155, 155, 156, 157,
      158, 159, 160, 160, 161, 162, 163, 163, 164, 165, 166, 167, 167, 168,
      169, 170, 170, 171, 172, 173, 173, 174, 175, 176, 176, 177, 178, 178,
      179, 180, 181, 181, 182, 183, 183, 184, 185, 185, 186, 187, 187, 188,
      189, 189, 190, 191, 192, 192, 193, 193, 194, 195, 195, 196, 197, 197,
      198, 199, 199, 200, 201, 201, 202, 203, 203, 204, 204, 205, 206, 206,
      207, 208, 208, 209, 209, 210, 211, 211, 212, 212, 213, 214, 214, 215,
      215, 216, 217, 217, 218, 218, 219, 219, 220, 221, 221, 222, 222, 223,
      224, 224, 225, 225, 226, 226, 227, 227, 228, 229, 229, 230, 230, 231,
      231, 232, 232, 233, 234, 234, 235, 235, 236, 236, 237, 237, 238, 238,
      239, 240, 240, 241, 241, 242, 242, 243, 243, 244, 244, 245, 245, 246,
      246, 247, 247, 248, 248, 249, 249, 250, 250, 251, 251, 252, 252, 253,
      253, 254, 254, 255
    };

    uint32_t xn;

    if (x >= 0x10000)
      if (x >= 0x1000000)
	if (x >= 0x10000000)
	  if (x >= 0x40000000) {
	    if (x >= 65535UL*65535UL)
	      return 65535;
	    xn = sqq_table[x>>24] << 8;
	  } else
	    xn = sqq_table[x>>22] << 7;
	else
	  if (x >= 0x4000000)
	    xn = sqq_table[x>>20] << 6;
	  else
	    xn = sqq_table[x>>18] << 5;
      else {
	if (x >= 0x100000)
	  if (x >= 0x400000)
	    xn = sqq_table[x>>16] << 4;
	  else
	    xn = sqq_table[x>>14] << 3;
	else
	  if (x >= 0x40000)
	    xn = sqq_table[x>>12] << 2;
	  else
	    xn = sqq_table[x>>10] << 1;
	
	goto nr1;
      }
    else
      if (x >= 0x100) {
	if (x >= 0x1000)
	  if (x >= 0x4000)
	    xn = (sqq_table[x>>8] >> 0) + 1;
	  else
	    xn = (sqq_table[x>>6] >> 1) + 1;
	else
	  if (x >= 0x400)
	    xn = (sqq_table[x>>4] >> 2) + 1;
	  else
	    xn = (sqq_table[x>>2] >> 3) + 1;
	
	goto adj;
      } else
	return sqq_table[x] >> 4;
    
    xn = (xn + 1 + x / xn) / 2;

  nr1:
    xn = (xn + 1 + x / xn) / 2;

  adj:
    if (xn * xn > x)
      --xn;
    
    return xn;
  }
  
} // namespace sl

namespace sl {

  static short ABDPY_shortdouble(double x) {
    static union {double number; short adr[sizeof(double)/sizeof(short)]; } val;     /* for lazy evaluation*/
    val.number=x;
    return val.adr[0];
  }

  
  const double ABDPY_2EXP53=9007199254740992.0;     /* 2^53*/

  static int ABDPY__not_lazy_det2x2(double x1, double y1, double x2, double y2) {
    int sign;
    double swap;
    double k;					/* quotient*/
    sign=1;
    /* testing null entries*/
    if ((x1 == 0.0)||(y2 == 0.0)) {                         /* det = -x2 y1*/
      if ((y1 == 0.0)||(x2 == 0.0)) {                     /* det = 0.0*/
        return 0;
      } else if (y1>0) {
        if (x2>0) {
          return -sign;
        } else {
          return  sign;
	}
      } else { 
        if (x2>0) {
          return  sign;
	} else {
          return -sign;
	}
      }
    }
    if ((y1 == 0.0)||(x2 == 0.0)) {                         /* det =  x1 y2*/
      if (y2>0) {
        if (x1>0) {
          return  sign;
        } else {
          return -sign;
	}
      } else { 
        if (x1>0) {
          return -sign;
        } else {
          return  sign;
	}
      }
    }

    /* making y coordinates positive and permuting the entries*/
    /* so that y2 is the biggest one*/
    {
      if ( 0.0 < y1 ) {                                /* 0.0 < y1*/
        if ( 0.0 < y2 ) {                               /* 0.0 < y1, y2*/
          if ( y1 <= y2 ){                               /* 0.0 < y1 <= y2*/
            /* everything ok*/
            ;
          }else{                                         /* 0.0 < y2 < y1*/
            /* permuting U1 U2*/
            sign = -sign;
            swap=x1; x1=x2; x2=swap;
            swap=y1; y1=y2; y2=swap;
          } } else {                                       /* y2 < 0.0 < y1*/
            if ( y1 <= -y2 ){                              /* 0.0 < y1 <= -y2*/
              /* changing U2 sign*/
              sign = -sign;
              x2 = -x2; y2 = -y2;
            }else{                                         /* 0.0 < -y2 < y1*/
              /* changing U2 sign and permuting U1 U2*/
              swap=x1; x1=-x2; x2=swap;
              swap=y1; y1=-y2; y2=swap;
            }} } else {                                       /* 0.0 < -y1*/
              if ( 0.0 < y2 ) {                               /* 0.0 < -y1, y2*/
                if ( -y1 <= y2 ){                              /* 0.0 < -y1 <= y2*/
                  /* changing U1 sign*/
                  sign = -sign;
                  x1 = -x1; y1 = -y1;
                }else{                                         /* 0.0 < y2 < -y1*/
                  /* changing U1 sign and permuting U1 U2*/
                  swap=-x1; x1=x2; x2=swap;
                  swap=-y1; y1=y2; y2=swap;
                } } else {                                       /* y2 < 0.0 < -y1*/
                  if ( y1 >= y2 ){                               /* 0.0 < -y1 <= -y2*/
                    /* changing U1 U2 sign*/
                    x1 = -x1; y1 = -y1;
                    x2 = -x2; y2 = -y2;
                    ;
                  }else{                                         /* 0.0 < -y2 < -y1*/
                    /* changing U1 U2 sign and permuting U1 U2*/
                    sign = -sign;
                    swap=-x1; x1=-x2; x2=swap;
                    swap=-y1; y1=-y2; y2=swap;
                  }} }
      /* making x coordinates positive*/
      /* if |x2| < |x1| one can conclude*/
      if ( 0.0 < x1 ) {                                /* 0.0 < x1*/
        if ( 0.0 < x2 ) {                               /* 0.0 < x1, x2*/
          if ( x1 <= x2 ){                               /* 0.0 < x1 < x2*/
            /* everything ok*/
            ;
          }else{                                         /* 0.0 < x2 < x1*/
            return sign;
          } } else {                                       /* 0.0 < x1, -x2*/
            return sign;
          }} else {                                         /* 0.0 < -x1*/
            if ( 0.0 < x2 ) {                               /* 0.0 < -x1, x2*/
              return -sign;
            } else {                                         /* x2 < 0.0 < -x1*/
              if ( x1 >= x2 ){                               /* 0.0 < -x1 <= -x2*/
                /* changing x1 x2 sign*/
                sign = -sign;
                x1 = -x1; x2 = -x2;
                ;
              }else{                                         /* 0.0 < -x2 < -x1*/
                return -sign;
              }} }
    }	
    /* all entries strictly positive   x1 <= x2 and y1 <= y2*/
    while (1)
      {
        k = floor(x2/x1); x2 -=k*x1; y2 -=k*y1;
        /* testing if R (new U2) is in U1 rectangle*/
        if (y2<0.0) return -sign;
        if (y2>y1 ) return  sign;
        /* finding R'*/
        if ( x1 > x2+x2 ) {                                         /* 0 is candidate*/
          if ( y1 < y2 + y2  ) {
            return sign ;
	  }
        } else {                                                      /* u1 is candidate*/
          if ( y1 > y2 + y2  ) {
            return -sign;
          } else {
	    x2=x1-x2;y2=y1-y2;sign=-sign; /* R' = u1-R*/
	  }
	}
   
        if (y2==0.0) return (x2==0.0) ? 0 : -sign;
        if (x2==0.0) return  sign;
        /* exchange 1 and 2 role.*/
        k = floor(x1/x2); x1 -=k*x2; y1 -=k*y2;
        /* testing if R (new U1) is in U2 rectangle*/
        if (y1<0.0) return  sign;
        if (y1>y2 ) return -sign;
        /* finding R'*/
        if ( x2 > x1+x1 ) {                                        /* 0 is candidate*/
          if ( y2 < y1 + y1  ) {
            return -sign ;
          } 
        } else {                                                      /* u2 is candidate*/
          if ( y2 > y1 + y1  ) {
            return  sign;
          } else {
	    x1=x2-x1;y1=y2-y1;sign=-sign;  /* R' = u2-R*/
	  }
	}
        if (y1==0.0) return (x1==0.0) ? 0 :  sign;
        if (x1==0.0) return -sign;
      }
  }

  int ABDPY__det2x2(double& x1, double& y1, double& x2, double& y2)
  {
    int sign;
    double swap;
    double k;					/* quotient*/
    sign=1;
    /* testing null entries*/
    if ((x1 == 0.0)||(y2 == 0.0)) {                          /* det = -x2 y1*/
      if ((y1 == 0.0)||(x2 == 0.0)) {                     /* det = 0.0*/
        return 0;
      } else if (y1>0) {
        if (x2>0) {
          return -sign;
        } else {
          return  sign;
	}
      } else { 
        if (x2>0) {
          return  sign;
	} else {
          return -sign;
	}
      }
    }
    if ((y1 == 0.0)||(x2 == 0.0)) {                         /* det =  x1 y2*/
      if (y2>0) {
        if (x1>0) {
          return  sign;
        } else {
          return -sign;
	}
      } else { 
        if (x1>0) {
          return -sign;
        } else {
          return  sign;
	}
      }
    }
    /* making y coordinates positive and permuting the entries*/
    /* so that y2 is the biggest one*/
    {
      if ( 0.0 < y1 ) {                                /* 0.0 < y1*/
        if ( 0.0 < y2 ) {                               /* 0.0 < y1, y2*/
          if ( y1 <= y2 ){                               /* 0.0 < y1 <= y2*/
            /* everything ok*/
            ;
          }else{                                         /* 0.0 < y2 < y1*/
            /* permuting U1 U2*/
            sign = -sign;
            swap=x1; x1=x2; x2=swap;
            swap=y1; y1=y2; y2=swap;
          } } else {                                       /* y2 < 0.0 < y1*/
            if ( y1 <= -y2 ){                              /* 0.0 < y1 <= -y2*/
              /* changing U2 sign*/
              sign = -sign;
              x2 = -x2; y2 = -y2;
            }else{                                         /* 0.0 < -y2 < y1*/
              /* changing U2 sign and permuting U1 U2*/
              swap=x1; x1=-x2; x2=swap;
              swap=y1; y1=-y2; y2=swap;
            }} } else {                                       /* 0.0 < -y1*/
              if ( 0.0 < y2 ) {                               /* 0.0 < -y1, y2*/
                if ( -y1 <= y2 ){                              /* 0.0 < -y1 <= y2*/
                  /* changing U1 sign*/
                  sign = -sign;
                  x1 = -x1; y1 = -y1;
                }else{                                         /* 0.0 < y2 < -y1*/
                  /* changing U1 sign and permuting U1 U2*/
                  swap=-x1; x1=x2; x2=swap;
                  swap=-y1; y1=y2; y2=swap;
                } } else {                                       /* y2 < 0.0 < -y1*/
                  if ( y1 >= y2 ){                               /* 0.0 < -y1 <= -y2*/
                    /* changing U1 U2 sign*/
                    x1 = -x1; y1 = -y1;
                    x2 = -x2; y2 = -y2;
                    ;
                  }else{                                         /* 0.0 < -y2 < -y1*/
                    /* changing U1 U2 sign and permuting U1 U2*/
                    sign = -sign;
                    swap=-x1; x1=-x2; x2=swap;
                    swap=-y1; y1=-y2; y2=swap;
                  }} }
      /* making x coordinates positive*/
      /* if |x2| < |x1| one can conclude*/
      if ( 0.0 < x1 ) {                                /* 0.0 < x1*/
        if ( 0.0 < x2 ) {                               /* 0.0 < x1, x2*/
          if ( x1 <= x2 ){                               /* 0.0 < x1 < x2*/
            /* everything ok*/
            ;
          }else{                                         /* 0.0 < x2 < x1*/
            return sign;
          } } else {                                       /* 0.0 < x1, -x2*/
            return sign;
          }} else {                                         /* 0.0 < -x1*/
            if ( 0.0 < x2 ) {                               /* 0.0 < -x1, x2*/
              return -sign;
            } else {                                         /* x2 < 0.0 < -x1*/
              if ( x1 >= x2 ){                               /* 0.0 < -x1 <= -x2*/
                /* changing x1 x2 sign*/
                sign = -sign;
                x1 = -x1; x2 = -x2;
                ;
              }else{                                         /* 0.0 < -x2 < -x1*/
                return -sign;
              }} }
    }	
    /* all entries strictly positive   x1 <= x2 and y1 <= y2*/
    while (1)
      {
        k = floor(x2/x1); x2 -=k*x1; y2 -=k*y1;
        /* testing if R (new U2) is in U1 rectangle*/
        if (y2<0.0) return -sign;
        if (y2>y1 ) return  sign;
        /* finding R'*/
        if ( x1 > x2+x2 ) {                                        /* 0 is candidate*/
          if ( y1 < y2 + y2  ) {
            return sign ;
          }
        } else {                                                      /* u1 is candidate*/
          if ( y1 > y2 + y2  ) {
            return -sign;
          } else {
	    x2=x1-x2;y2=y1-y2;sign=-sign;  /* R' = u1-R*/
	  }
	}
        if (y2==0.0) return (x2==0.0) ? 0 : -sign;
        if (x2==0.0) return  sign;
        /* exchange 1 and 2 role.*/
        k = floor(x1/x2); x1 -=k*x2; y1 -=k*y2;
        /* testing if R (new U1) is in U2 rectangle*/
        if (y1<0.0) return  sign;
        if (y1>y2 ) return -sign;
        /* finding R'*/
        if ( x2 > x1+x1 ) {                                        /* 0 is candidate*/
          if ( y2 < y1 + y1  ) {
            return -sign ;
          }
        } else {                                                      /* u2 is candidate*/
          if ( y2 > y1 + y1  ) {
            return  sign;
          } else {                                                   
	    x1=x2-x1;y1=y2-y1;sign=-sign; /* R' = u2-R*/ 
	  }
	}
        if (y1==0.0) return (x1==0.0) ? 0 :  sign;
        if (x1==0.0) return -sign;
        /* LAZY evaluation*/
        swap = x1*y2; k = x2*y1;
        if (swap > k) return  sign;
        if (swap < k) return -sign;
      }
  }








  int ABDPY_det2x2(double x1, double y1, double x2, double y2)
  {
    double p,q;
    /* LAZY EVALUATION*/
    p =x1*y2;
    q =x2*y1;
    if (p > q) return 1;
    if (p < q) return -1;
    return ABDPY__det2x2(x1,y1,x2,y2);
  }
  
  int ABDPY__not_lazy_det3x3(double x1, double y1, double z1, 
                             double x2, double y2, double z2,
                             double x3, double y3, double z3)
  {
    register double swap;
    static   double xx,yy,zz,xxx,yyy/*,zzz*/;
    static   int sign;            /* difference between the sign of the original*/
    /* determinant and the current one (value 1 or -1)*/
    static   int m1,m2,m3;        /* sign of minors*/
    static   double k;            /* quotient in some special cases (null minors)*/
    static   double xv1,yv1,zv1;  /* v1 vector*/
    static   double xv2,yv2,zv2;  /* v2 vector*/
    static   int u_to_v;    	  /* difference between ui and vi*/
    static   double cccx[106];    /* stack of C(l) points*/
    static   double cccy[106];    /* stack of C(l) points*/
    static   double cccz[106];    /* stack of C(l) points*/
    register double* cx;		  /* current place in the stack*/
    register double* cy;		  /* current place in the stack*/
    register double* cz;		  /* current place in the stack*/
    static   double ccx,ccy,ccz;  /* C(l) in step 1.4*/

    sign=1;
    /* first iteration*/
    {
      /* making z coordinates positive and permuting the entries*/
      /* so that z3 is the biggest one*/
      {
        if ( 0.0 <= z1 ){                                   /* 0 <= z1*/
          if ( 0.0 <= z2 ){                                  /* 0 <= z1, z2*/
            if ( z1 <= z3 ){                                /* 0 <= z1, z2  and z1<=z3*/
              if ( z2 <= z3 ){                               /* 0 <= z1, z2 <= z3*/
                /* everything ok*/
                ;
              }else{                                         /* 0 <= z1 <= z3 < = z2*/
                /* permuting U2 U3*/
                sign = -sign;
                swap=x2; x2=x3; x3=swap;
                swap=y2; y2=y3; y3=swap;
                swap=z2; z2=z3; z3=swap;
              }}else{                                         /* 0 <= z1, z2  and z3<z1*/
                if ( 0.0 <= z3 ) {                               /* 0 <= z2, z3  and z3<z1*/
                  if ( z1 <= z2 ){                              /* 0 <= z3 < z1 <= z2*/
                    /* permuting U2 U3*/
                    sign = -sign;
                    swap=x2; x2=x3; x3=swap;
                    swap=y2; y2=y3; y3=swap;
                    swap=z2; z2=z3; z3=swap;
                  }else{                                        /* 0 <= z2, z3 < z1*/
                    /* permuting U1 U3*/
                    sign = -sign;
                    swap=x1; x1=x3; x3=swap;
                    swap=y1; y1=y3; y3=swap;
                    swap=z1; z1=z3; z3=swap;
                  }}else{                                        /* z3 < 0 <= z1, z2*/
                    if ( z1 <= -z3 ){                             /* 0 <= z1, z2 and z1 <= -z3*/
                      if ( z2 <= -z3 ){                             /* 0 <= z1, z2 <= -z3*/
                        /* changing U3 sign*/
                        sign = -sign;
                        x3 = -x3; y3 = -y3; z3 = -z3;
                      }else{                                        /* 0 <= z1 <= -z3 < z2*/
                        /* changing U3 sign and permuting U2 U3 */
                        swap=x2; x2=-x3; x3=swap;
                        swap=y2; y2=-y3; y3=swap;
                        swap=z2; z2=-z3; z3=swap;
                      }}else{                                        /* 0 <= z1, z2 and -z3 < z1*/
                        if ( z2 <= z1 ){                              /* 0 <= -z3, z2 <= z1*/
                          /* changing U3 sign and permuting U1 U3*/
                          swap=x1; x1=-x3; x3=swap;
                          swap=y1; y1=-y3; y3=swap;
                          swap=z1; z1=-z3; z3=swap;
                        }else{                                        /* 0 <= -z3 < z1 < z2*/
                          /* changing U3 sign and permuting U2 U3*/
                          swap=x2; x2=-x3; x3=swap;
                          swap=y2; y2=-y3; y3=swap;
                          swap=z2; z2=-z3; z3=swap;
                        }}}}}else{                                        /* z2 < 0 <= z1*/
                          if ( z1 <= z3 ){                                 /* z2 < 0 <= z1 <= z3*/
                            if ( -z2 <= z3 ){                               /* 0 <= z1, -z2 <= z3*/
                              /* changing U2 sign*/
                              sign = -sign;
                              x2 = -x2; y2 = -y2; z2 = -z2;
                            }else{                                          /* 0 <= z1 <= z3 < -z2*/
                              /* changing U2 sign and permuting U2 U3 */
                              swap=-x2; x2=x3; x3=swap;
                              swap=-y2; y2=y3; y3=swap;
                              swap=-z2; z2=z3; z3=swap;
                            }}else{                                           /* z2 < 0 <= z1 and z3<z1*/
                              if ( 0.0 <= z3 ) {                                 /* 0<=z1, -z2, z3 and z3<z1*/
                                if ( -z2 <= z1 ){                                /* 0 <= -z2, z3 <= z1*/
                                  /* changing U2 sign and permuting U1 U3 */
                                  x2 = -x2; y2 = -y2; z2 = -z2;
                                  swap=x1; x1=x3; x3=swap;
                                  swap=y1; y1=y3; y3=swap;
                                  swap=z1; z1=z3; z3=swap;
                                }else{                                           /* 0 <= z3 <= z1 < -z2*/
                                  /* changing U2 sign and permuting U2 U3 */
                                  swap=-x2; x2=x3; x3=swap;
                                  swap=-y2; y2=y3; y3=swap;
                                  swap=-z2; z2=z3; z3=swap;
                                }}else{                                           /* 0<=z1, -z2, -z3*/
                                  if ( z1 <= -z3) {                               /* z2 < 0 <= z1 <= -z3*/
                                    if ( z3 <= z2){                                /* 0 <= z1, -z2 <= -z3*/
                                      /* changing U2 U3 sign*/
                                      x2 = -x2; y2 = -y2; z2 = -z2;
                                      x3 = -x3; y3 = -y3; z3 = -z3;
                                    }else{                                         /* 0 <= z1 <= -z3 < -z2*/
                                      /* changing U2 U3 sign and permuting U2 U3 */
                                      sign = -sign;
                                      swap=-x2; x2=-x3; x3=swap;
                                      swap=-y2; y2=-y3; y3=swap;
                                      swap=-z2; z2=-z3; z3=swap;
                                    }}else {                                        /* z2 < 0 <= -z3 < z1*/
                                      if ( z1 <= -z2 ) {                             /* 0 <= -z3 < z1 <= -z2*/
					/* changing U2 U3 sign and permuting U2 U3 */
					sign = -sign;
					swap=-x2; x2=-x3; x3=swap;
					swap=-y2; y2=-y3; y3=swap;
					swap=-z2; z2=-z3; z3=swap;
                                      }else{                                         /* 0 <= -z3, -z2 < z1*/
					/* changing U2 U3 sign and permuting U1 U3 */
					sign = -sign;
					x2 = -x2; y2 = -y2; z2 = -z2;
					swap=x1; x1=-x3; x3=swap;
					swap=y1; y1=-y3; y3=swap;
					swap=z1; z1=-z3; z3=swap;
                                      }}}}}}else{                                          /* z1 < 0*/
                                        if ( 0.0 <= z2 ){                                  /* 0 <= -z1, z2*/
                                          if ( -z1 <= z3 ){                                /* 0 <= -z1, z2  and -z1<=z3*/
                                            if ( z2 <= z3 ){                               /* 0 <= -z1, z2 <= z3*/
                                              /* changing U1 sign*/
                                              sign = -sign;
                                              x1 = -x1; y1 = -y1; z1 = -z1;
                                            }else{                                         /* 0 <= -z1 <= z3 < = z2*/
                                              /* changing U1 sign permuting U2 U3*/
                                              x1 = -x1; y1 = -y1; z1 = -z1;
                                              swap=x2; x2=x3; x3=swap;
                                              swap=y2; y2=y3; y3=swap;
                                              swap=z2; z2=z3; z3=swap;
                                            }}else{                                         /* 0 <= -z1, z2  and z3<-z1*/
                                              if ( 0.0 <= z3 ) {                               /* 0 <= z2, z3  and z3<-z1*/
                                                if ( -z1 <= z2 ){                              /* 0 <= z3 < -z1 <= z2*/
                                                  /* changing U1 sign permuting U2 U3*/
                                                  x1 = -x1; y1 = -y1; z1 = -z1;
                                                  swap=x2; x2=x3; x3=swap;
                                                  swap=y2; y2=y3; y3=swap;
                                                  swap=z2; z2=z3; z3=swap;
                                                }else{                                        /* 0 <= z2, z3 < -z1*/
                                                  /* changing U1 sign permuting U1 U3*/
                                                  swap=-x1; x1=x3; x3=swap;
                                                  swap=-y1; y1=y3; y3=swap;
                                                  swap=-z1; z1=z3; z3=swap;
                                                }}else{                                        /* z3 < 0 <= -z1, z2*/
                                                  if ( -z1 <= -z3 ){                             /* 0 <=-z1,z2 and -z1 <= -z3*/
                                                    if ( z2 <= -z3 ){                             /* 0 <= -z1, z2 <= -z3*/
                                                      /* changing U1 U3 sign*/
                                                      x1 = -x1; y1 = -y1; z1 = -z1;
                                                      x3 = -x3; y3 = -y3; z3 = -z3;
                                                    }else{                                        /* 0 <= -z1 <= -z3 < z2*/
                                                      /* changing U1 U3 sign and permuting U2 U3 */
                                                      sign = -sign;
                                                      x1 = -x1; y1 = -y1; z1 = -z1;
                                                      swap=x2; x2=-x3; x3=swap;
                                                      swap=y2; y2=-y3; y3=swap;
                                                      swap=z2; z2=-z3; z3=swap;
                                                    }}else{                                        /* 0 <=-z1,z2 and -z3 < -z1*/
                                                      if ( z2 <= -z1 ){                              /* 0 <= -z3, z2 <= -z1*/
                                                        /* changing U1 U3 sign and permuting U1 U3*/
                                                        sign = -sign;
                                                        swap=-x1; x1=-x3; x3=swap;
                                                        swap=-y1; y1=-y3; y3=swap;
                                                        swap=-z1; z1=-z3; z3=swap;
                                                      }else{                                        /* 0 <= -z3 < -z1 < z2*/
                                                        /* changing U1 U3 sign and permuting U2 U3*/
                                                        sign = -sign;
                                                        x1 = -x1; y1 = -y1; z1 = -z1;
                                                        swap=x2; x2=-x3; x3=swap;
                                                        swap=y2; y2=-y3; y3=swap;
                                                        swap=z2; z2=-z3; z3=swap;
                                                      }}}}}else{                                        /* z2 < 0 <= -z1*/
                                                        if ( -z1 <= z3 ){                                 /* z2 < 0 <= -z1 <= z3*/
                                                          if ( -z2 <= z3 ){                               /* 0 <= -z1, -z2 <= z3*/
                                                            /* changing U1U2 sign*/
                                                            x1 = -x1; y1 = -y1; z1 = -z1;
                                                            x2 = -x2; y2 = -y2; z2 = -z2;
                                                          }else{                                          /* 0 <= -z1 <= z3 < -z2*/
                                                            /* changing U1 U2 sign and permuting U2 U3 */
                                                            sign = -sign;
                                                            x1 = -x1; y1 = -y1; z1 = -z1;
                                                            swap=-x2; x2=x3; x3=swap;
                                                            swap=-y2; y2=y3; y3=swap;
                                                            swap=-z2; z2=z3; z3=swap;
                                                          }}else{                                           /* z2<0 <= -z1 and z3<-z1*/
                                                            if ( 0.0 <= z3 ) {                                 /* 0<=-z1,-z2,z3 and z3<-z1*/
                                                              if ( z2 >= z1 ){                                 /* 0 <= -z2, z3 <= -z1*/
                                                                /* changing U1 U2 sign and permuting U1 U3 */
                                                                sign = -sign;
                                                                x2 = -x2; y2 = -y2; z2 = -z2;
                                                                swap=-x1; x1=x3; x3=swap;
                                                                swap=-y1; y1=y3; y3=swap;
                                                                swap=-z1; z1=z3; z3=swap;
                                                              }else{                                           /* 0 <= z3 <= -z1 < -z2*/
                                                                /* changing U1 U2 sign and permuting U2 U3 */
                                                                sign = -sign;
                                                                x1 = -x1; y1 = -y1; z1 = -z1;
                                                                swap=-x2; x2=x3; x3=swap;
                                                                swap=-y2; y2=y3; y3=swap;
                                                                swap=-z2; z2=z3; z3=swap;
                                                              }}else{                                           /* 0<=-z1, -z2, -z3*/
                                                                if ( -z1 <= -z3) {                               /* z2 < 0 <= -z1 <= -z3*/
                                                                  if ( z3 <= z2){                                /* 0 <= -z1, -z2 <= -z3*/
                                                                    /* changing U1 U2 U3 sign*/
                                                                    sign = -sign;
                                                                    x1 = -x1; y1 = -y1; z1 = -z1;
                                                                    x2 = -x2; y2 = -y2; z2 = -z2;
                                                                    x3 = -x3; y3 = -y3; z3 = -z3;
                                                                  }else{                                         /* 0 <= -z1 <= -z3 < -z2*/
                                                                    /* changing U1 U2 U3 sign and permuting U2 U3 */
                                                                    x1 = -x1; y1 = -y1; z1 = -z1;
                                                                    swap=-x2; x2=-x3; x3=swap;
                                                                    swap=-y2; y2=-y3; y3=swap;
                                                                    swap=-z2; z2=-z3; z3=swap;
                                                                  }}else {                                        /* z2 < 0 <= -z3 < -z1*/
                                                                    if ( -z1 <= -z2 ) {                             /* 0 <= -z3 < -z1 <= -z2*/
                                                                      /* changing U1 U2 U3 sign and permuting U2 U3 */
                                                                      x1 = -x1; y1 = -y1; z1 = -z1;
                                                                      swap=-x2; x2=-x3; x3=swap;
                                                                      swap=-y2; y2=-y3; y3=swap;
                                                                      swap=-z2; z2=-z3; z3=swap;
                                                                    }else{                                         /* 0 <= -z3, -z2 < -z1*/
                                                                      /* changing U1 U2 U3 sign and permuting U1 U3 */
                                                                      x2 = -x2; y2 = -y2; z2 = -z2;
                                                                      swap=-x1; x1=-x3; x3=swap;
                                                                      swap=-y1; y1=-y3; y3=swap;
                                                                      swap=-z1; z1=-z3; z3=swap;
                                                                    }}}}}}}
      if ((m3=ABDPY_det2x2(x1,y1,x2,y2))==0)
        {                  /* the lattice is undefined, but z3 can be replaced by 0.0*/
          z3=0.0;
          /* redo the permutation*/
          if ( z1 <= z3 ){                                /* 0 <= z1, z2  and z1<=z3*/
            if ( z2 <= z3 ){                               /* 0 <= z1, z2 <= z3*/
              /* everything ok cannot arrive!*/
              ;
            }else{                                         /* 0 <= z1 <= z3 < = z2*/
              /* permuting U2 U3*/
              sign = -sign;
              swap=x2; x2=x3; x3=swap;
              swap=y2; y2=y3; y3=swap;
              swap=z2; z2=z3; z3=swap;
            }}else{                                         /* 0 <= z1, z2  and z3<z1*/
              if ( z1 <= z2 ){                              /* 0 <= z3 < z1 <= z2*/
                /* permuting U2 U3*/
                sign = -sign;
                swap=x2; x2=x3; x3=swap;
                swap=y2; y2=y3; y3=swap;
                swap=z2; z2=z3; z3=swap;
              }else{                                        /* 0 <= z2, z3 < z1*/
                /* permuting U1 U3*/
                sign = -sign;
                swap=x1; x1=x3; x3=swap;
                swap=y1; y1=y3; y3=swap;
                swap=z1; z1=z3; z3=swap;
              }
            }
          if ((m3=ABDPY_det2x2(x1,y1,x2,y2))==0)
            {                 /* the lattice is twice undifined, z3 can be replaced by 0.0*/
              z3=0.0;
              /* now at least two z values are null*/
              if (sign == 1)
                if ( z1 == 0.0 )									/* z1 = z3 = 0 <= z2*/
                  if ( z2 == 0.0)
                    return 0;
                  else
                    return   ABDPY_det2x2(x3,y3,x1,y1) ;
                else								              	/* z2 = z3 = 0 < z1*/
                  return   ABDPY_det2x2(x2,y2,x3,y3) ;
              else
                if ( z1 == 0.0 )									/* z1 = z3 = 0 <= z2*/
                  if ( z2 == 0.0)
                    return 0;
                  else
                    return   ABDPY_det2x2(x1,y1,x3,y3) ;
                else								              	/* z2 = z3 = 0 < z1*/
                  return   ABDPY_det2x2(x3,y3,x2,y2) ;
            }
        }
    }	
    while(1)
      {
        if (z3 == 0.0) return 0;                            /* z1=z2=z3=0.0*/
        /* compute the remaining minors*/
        m1 = ABDPY_det2x2(x2,y2,x3,y3);
        m2 = ABDPY_det2x2(x3,y3,x1,y1);
        /* if the three minors have the same sign*/
        if ( ( m1 >= 0 ) && ( m2 >= 0) && (m3 >= 0)) return sign ;
        if ( ( m1 <= 0 ) && ( m2 <= 0) && (m3 <= 0)) return -sign ;
        /* if one minor is null the corresponding z value can be replaced by 0.0*/
        /* m3==0 already tested*/
        if ( m1 == 0 )
          {
            z1 = 0.0;
            if (x2 != 0.0)
              {
		k = floor(x3/x2);    /* as in 2D case k can be computed directly*/
		x3 -=k*x2;
		y3 -=k*y2;
		z3 -=k*z2;
              }
            else       /* x2=0.0*/
              {
		k = floor(y3/y2);    /* as in 2D case k can be computed directly*/
		y3 -=k*y2;
		z3 -=k*z2;
              }
          }
        else if ( m2 == 0 )            /* m1 != 0*/
          {
            z2 = 0.0;
            if (x1 != 0.0)
              {
		k = floor(x3/x1);    /* as in 2D case k can be computed directly*/
		x3 -=k*x1;
		y3 -=k*y1;
		z3 -=k*z1;
              }
            else       /* x1=0.0*/
              {
		k = floor(y3/y1);    /* as in 2D case k can be computed directly*/
		y3 -=k*y1;
		z3 -=k*z1;
              }
          }
        else                           /* m1 != 0 and m2 != 0 and m3 != 0*/
          {													/* step 1 : compute R*/
            /* substep 1.1*/
            /* computing v1 v2 such that v1 + v2 and u3 is in the positive wedge v1 v2*/
            if ( m3 > 0 )                   /*           m3>0*/
              if ( m1 > 0 )                  /* m1>0      m3>0*/
                if ( m2 > 0 )                 /* m1>0 m2>0 m3>0*/
                  {
                    u_to_v = 1;
                    xv1 = -x1; yv1 = -y1; zv1 = -z1;
                    xv2 = -x2; yv2 = -y2; zv2 = -z2;
                  }
                else                          /* m1>0 m2<0 m3>0*/
                  {
                    u_to_v = 2;
                    xv1 =  x2; yv1 =  y2; zv1 =  z2;
                    xv2 = -x1; yv2 = -y1; zv2 = -z1;
                  }
              else                           /* m1<0      m3>0*/
                if ( m2 > 0 )                 /* m1<0 m2>0 m3>0*/
                  {
                    u_to_v = 3;
                    xv1 = -x2; yv1 = -y2; zv1 = -z2;
                    xv2 =  x1; yv2 =  y1; zv2 =  z1;
                  }
                else                          /* m1<0 m2<0 m3>0*/
                  {
                    u_to_v = 4;
                    xv1 =  x1; yv1 =  y1; zv1 =  z1;
                    xv2 =  x2; yv2 =  y2; zv2 =  z2;
                  }
            else                            /*           m3<0*/
              if ( m1 > 0 )                  /* m1>0      m3<0*/
                if ( m2 > 0 )                 /* m1>0 m2>0 m3<0*/
                  {
                    u_to_v = 5;
                    xv1 =  x2; yv1 =  y2; zv1 =  z2;
                    xv2 =  x1; yv2 =  y1; zv2 =  z1;
                  }
                else                          /* m1>0 m2<0 m3<0*/
                  {
                    u_to_v = 6;
                    xv1 =  x1; yv1 =  y1; zv1 =  z1;
                    xv2 = -x2; yv2 = -y2; zv2 = -z2;
                  }
              else                           /* m1<0      m3<0*/
                if ( m2 > 0 )                 /* m1<0 m2>0 m3<0*/
                  {
                    u_to_v = 7;
                    xv1 = -x1; yv1 = -y1; zv1 = -z1;
                    xv2 =  x2; yv2 =  y2; zv2 =  z2;
                  }
                else                          /* m1<0 m2<0 m3<0*/
                  {
                    u_to_v = 8;
                    xv1 = -x2; yv1 = -y2; zv1 = -z2;
                    xv2 = -x1; yv2 = -y1; zv2 = -z1;
                  }
            /* substep 1.2*/
            /* verify if u3 is in the base cell of lattice v1 v2*/
            if ((  ABDPY_det2x2(xv1, yv1, x3-xv2, y3-yv2) > 0 ) ||
                (  ABDPY_det2x2(x3-xv1, y3-yv1, xv2, yv2) > 0 ))
              {
                if (  ABDPY_det2x2(x3, y3, xv1+xv2, yv1+yv2) > 0 )
                  /* u3 in the wedge v1 v1+v2*/
                  {
                    /* substep 1.3*/
                    *(cx=cccx) = xv1;
                    *(cy=cccy) = yv1;
                    *(cz=cccz) = zv1;
                    while (1)
                      { 
                        zz=*cz+zv2;
                        if ((*cz<0)&&(zz<0)) return (m3>0.0) ? sign : -sign ;
                        if ((*cz>ABDPY_2EXP53)&&(zz>ABDPY_2EXP53))
                          return (m3>0.0) ? -sign : sign;
                        xx=*cx+*cx; 
                        yy=*cy+*cy; 
                        if ( ABDPY_det2x2(xv2, yv2, x3-xx, y3-yy) >= 0)
                          break; /* goto 1.4*/
                        xxx = xx + xv2;
                        yyy = yy + yv2;
                        if ( ABDPY_det2x2(x3, y3, xxx, yyy) >= 0)
                          { *++cx =  xx; *++cy =  yy; zz= *cz+*cz; *++cz = zz;}
                        else
                          { *++cx = xxx; *++cy = yyy; zz+= *cz;   *++cz = zz;}
                      }
                    /* substep 1.4*/
                    ccx = *cx;
                    ccy = *cy;
                    ccz = *cz;
                    while (--cx >= cccx)
                      { 
                        --cy;--cz;
                        zz=ccz+zv2;
                        if ((ccz<0)&&(zz<0)) return (m3>0.0) ? sign : -sign ;
                        if ((ccz>ABDPY_2EXP53)&&(zz>ABDPY_2EXP53))
                          return (m3>0.0) ? -sign : sign;
                        xx=ccx+*cx; 
                        yy=ccy+*cy; 
                        if ( ABDPY_det2x2(xv2, yv2, x3-xx, y3-yy) >= 0)
                          continue;
                        xxx = xx + xv2;
                        yyy = yy + yv2;
                        if ( ABDPY_det2x2(x3, y3, xxx, yyy) >= 0)
                          { ccx =  xx; ccy =  yy; ccz+= *cz;}
                        else
                          { ccx = xxx; ccy = yyy; ccz = zz+*cz;}
                      }
                    /* substep 1.5*/
                    xxx = ccx + xv2;
                    yyy = ccy + yv2;
                    if ( ABDPY_det2x2(xv1, yv1, x3-xxx, y3-yyy) > 0)
                      { ccx = xxx; ccy = yyy; ccz += zv2; }
                  }
                else
                  /* u3 in the wedge v1+v2 v2*/
                  {
                    /* substep 1.3*/
                    *(cx=cccx) = xv2;
                    *(cy=cccy) = yv2;
                    *(cz=cccz) = zv2;
                    while (1)
                      { 
                        zz=*cz+zv1;
                        if ((*cz<0)&&(zz<0)) return (m3>0.0) ? sign : -sign ;
                        if ((*cz>ABDPY_2EXP53)&&(zz>ABDPY_2EXP53))
                          return (m3>0.0) ? -sign : sign;
                        xx=*cx+*cx; 
                        yy=*cy+*cy; 
                        if ( ABDPY_det2x2(xv1, yv1, x3-xx, y3-yy) <= 0)
                          break; /* goto 1.4*/
                        xxx = xx + xv1;
                        yyy = yy + yv1;
                        if ( ABDPY_det2x2(x3, y3, xxx, yyy) <= 0)
                          { *++cx =  xx; *++cy =  yy; zz= *cz+*cz; *++cz = zz;}
                        else
                          { *++cx = xxx; *++cy = yyy; zz+= *cz; *++cz = zz;}
                      }
                    /* substep 1.4*/
                    ccx = *cx;
                    ccy = *cy;
                    ccz = *cz;
                    while (--cx >= cccx)
                      { 
                        --cy;--cz;
                        zz=ccz+zv1;
                        if ((ccz<0)&&(zz<0)) return (m3>0.0) ? sign : -sign ;
                        if ((ccz>ABDPY_2EXP53)&&(zz>ABDPY_2EXP53))
                          return (m3>0.0) ? -sign : sign;
                        xx=ccx+*cx; 
                        yy=ccy+*cy; 
                        if ( ABDPY_det2x2(xv1, yv1, x3-xx, y3-yy) <= 0)
                          continue;
                        xxx = xx + xv1;
                        yyy = yy + yv1;
                        if ( ABDPY_det2x2(x3, y3, xxx, yyy) <= 0)
                          { ccx =  xx; ccy =  yy; ccz+= *cz;}
                        else
                          { ccx = xxx; ccy = yyy; ccz = zz+*cz;}
                      }
                    /* substep 1.5*/
                    xxx = ccx + xv1;
                    yyy = ccy + yv1;
                    if ( ABDPY_det2x2(xv2, yv2, x3-xxx, y3-yyy) < 0)
                      { ccx = xxx; ccy = yyy; ccz += zv1; }
                  }	  
                switch (u_to_v)
                  {
                  case 1:
                  case 8:
                    ccx += xv1+xv2; ccy += yv1+yv2; ccz += zv1+zv2; break;
                  case 3:
                  case 7:
                    ccx += xv1; ccy += yv1; ccz += zv1; break;
                  case 2:
                  case 6:
                    ccx += xv2; ccy += yv2; ccz += zv2; break;
                  }
                x3 -= ccx; y3 -= ccy; z3 -= ccz;
              }
            else             /* u3 is originally in the base cell of lattice v1 v2*/
              switch (u_to_v)
                {
                case 1:
                case 8:
                  x3 -= xv1+xv2; y3 -= yv1+yv2; z3 -= zv1+zv2; break;
                case 3:
                case 7:
                  x3 -= xv1; y3 -= yv1; z3 -= zv1; break;
                case 2:
                case 6:
                  x3 -= xv2; y3 -= yv2; z3 -= zv2; break;
                }
          }
        /* substep 1.6*/
        if (z3<0.0) return (m3>0.0) ? -sign : sign;
        if (z3>z1+z2) return (m3>0.0) ? sign : -sign ;
        {
          /* step 2*/
          /* finding R' with z >= 0*/
          if ( ABDPY_det2x2(x2-x1, y2-y1, x3-x1, y3-y1) > 0 )
            {if ( ABDPY_det2x2(x3, y3, x1+x2, y1+y2) > 0 )
              {if ( m3 > 0.0 )                                              /* 0 or u1*/
                {if ( z1 > z2 )
                  if ( z3 > z1 ) return sign;
                  else {if ( z3 + z3 > z1 ) {sign=-sign;z3=z1-z3; x3=x1-x3; y3=y1-y3;}}
                else
                  if ( (ccz=z3 + z3) > z1 + z2 ) return sign;
                  else {if ( z3 > z1 )  {z3-=z1; x3-=x1; y3-=y1;}
                  else if ( ccz > z1 ) {sign=-sign;z3=z1-z3; x3=x1-x3; y3=y1-y3;}}
                }
              else                                                         /* u1+u2 or u2*/
                {if ( z2 < z1 )
                  if ( z3 < z2 ) return sign;
                  else {z3-=z2;x3-=x2;y3-=y2;
                  if ( z3+z3 > z1 ) {sign=-sign;z3=z1-z3; x3=x1-x3; y3=y1-y3;}}
                else
                  if ( z3 + z3 < z1 + z2 ) return sign;
                  else {if ( z3 < z2 ) {sign=-sign;z3=z2-z3; x3=x2-x3; y3=y2-y3;}
                  else{z3-=z2;x3-=x2;y3-=y2;
                  if ( z3+z3 > z1 ) {sign=-sign;z3=z1-z3; x3=x1-x3; y3=y1-y3;}}}
                }
              }
            else
              {if ( m3 > 0.0 )                                              /* 0 or u2*/
                {if ( z2 > z1 )
                  if ( z3 > z2 ) return sign;
                  else {if ( z3 + z3 > z2 ) {sign=-sign;z3=z2-z3; x3=x2-x3; y3=y2-y3;}}
                else
                  if ( (ccz=z3 + z3) > z1 + z2 ) return sign;
                  else {if ( z3 > z2 ) {z3-=z2;x3-=x2;y3-=y2;}
                  else if ( ccz > z2 ) {sign=-sign;z3=z2-z3; x3=x2-x3; y3=y2-y3;}}
                }
              else                                                         /* u1+u2 or u1*/
                {if ( z1 < z2 )
                  if ( z3 < z1 ) return sign;
                  else {z3-=z1;x3-=x1;y3-=y1;
                  if ( z3+z3 > z2 ) {sign=-sign;z3=z2-z3; x3=x2-x3; y3=y2-y3;}}
                else
                  if ( z3 + z3 < z1 + z2 ) return sign;
                  else {if ( z3 < z1 ) {sign=-sign;z3=z1-z3; x3=x1-x3; y3=y1-y3;}
                  else{z3-=z1;x3-=x1;y3-=y1;
                  if ( z3+z3 > z2 ) {sign=-sign;z3=z2-z3; x3=x2-x3; y3=y2-y3;}}}
                }
              }
            }
          else
            {if ( ABDPY_det2x2(x3, y3, x1+x2, y1+y2) > 0 )
              {if ( m3 > 0.0 )                                              /* u1+u2 or u1*/
                {if ( z1 < z2 )
                  if ( z3 < z1 ) return -sign;
                  else {z3-=z1;x3-=x1;y3-=y1;
                  if ( z3+z3 > z2 ) {sign=-sign;z3=z2-z3; x3=x2-x3; y3=y2-y3;}}
                else
                  if ( z3 + z3 < z1 + z2 ) return -sign;
                  else {if ( z3 < z1 ) {sign=-sign;z3=z1-z3; x3=x1-x3; y3=y1-y3;}
                  else{z3-=z1;x3-=x1;y3-=y1;
                  if ( z3+z3 > z2 ) {sign=-sign;z3=z2-z3; x3=x2-x3; y3=y2-y3;}}}
                }
              else                                                         /* 0 or u2*/
                {if ( z2 > z1 )
                  if ( z3 > z2 ) return -sign;
                  else {if ( z3 + z3 > z2 ) {sign=-sign;z3=z2-z3; x3=x2-x3; y3=y2-y3;}}
                else
                  if ( (ccz=z3 + z3) > z1 + z2 ) return -sign;
                  else {if ( z3 > z2 ) {z3-=z2;x3-=x2;y3-=y2;}
                  else if ( ccz > z2 ) {sign=-sign;z3=z2-z3; x3=x2-x3; y3=y2-y3;}}
                }
              }
            else
              {if ( m3 > 0.0 )                                              /* u1+u2 or u2*/
                {if ( z2 < z1 )
                  if ( z3 < z2 ) return -sign;
                  else {z3-=z2;x3-=x2;y3-=y2;
                  if ( z3+z3 > z1 ) {sign=-sign;z3=z1-z3; x3=x1-x3; y3=y1-y3;}}
                else
                  if ( z3 + z3 < z1 + z2 ) return -sign;
                  else {if ( z3 < z2 ) {sign=-sign;z3=z2-z3; x3=x2-x3; y3=y2-y3;}
                  else{z3-=z2;x3-=x2;y3-=y2;
                  if ( z3+z3 > z1 ) {sign=-sign;z3=z1-z3; x3=x1-x3; y3=y1-y3;}}}
                }
              else                                                         /* 0 or u1*/
                {if ( z1 > z2 )
                  if ( z3 > z1 ) return -sign;
                  else {if ( z3 + z3 > z1 ) {sign=-sign;z3=z1-z3; x3=x1-x3; y3=y1-y3;}}
                else
                  if ( (ccz=z3 + z3) > z1 + z2 ) return -sign;
                  else {if ( z3 > z1 )  {z3-=z1; x3-=x1; y3-=y1;}
                  else if ( ccz > z1 ) {sign=-sign;z3=z1-z3; x3=x1-x3; y3=y1-y3;}}
                }
              }
            }
        } /* end R'*/
        /* further iteration  0 <= z1,z2,z3 <= max(z1,z2)/2*/
        {
          /* making z coordinates positive and permuting the entries*/
          /* so that z3 is the biggest one*/
          if ( z1 <= z3 ){                                /* 0 <= z1, z2  and z1<=z3*/
            if ( z2 <= z3 ){                               /* 0 <= z1, z2 <= z3*/
              /* everything ok cannot arrive!*/
              ;
            }else{                                         /* 0 <= z1 <= z3 < = z2*/
              /* permuting U2 U3*/
              sign = -sign;
              swap=x2; x2=x3; x3=swap;
              swap=y2; y2=y3; y3=swap;
              swap=z2; z2=z3; z3=swap;
            }}else{                                         /* 0 <= z1, z2  and z3<z1*/
              if ( z1 <= z2 ){                              /* 0 <= z3 < z1 <= z2*/
                /* permuting U2 U3*/
                sign = -sign;
                swap=x2; x2=x3; x3=swap;
                swap=y2; y2=y3; y3=swap;
                swap=z2; z2=z3; z3=swap;
              }else{                                        /* 0 <= z2, z3 < z1*/
                /* permuting U1 U3*/
                sign = -sign;
                swap=x1; x1=x3; x3=swap;
                swap=y1; y1=y3; y3=swap;
                swap=z1; z1=z3; z3=swap;
              }
            }}
        if ((m3=ABDPY_det2x2(x1,y1,x2,y2))==0)
          {                  /* the lattice is undefined, but z3 can be replaced by 0.0*/
            z3=0.0;
            /* redo the permutation*/
            if ( z1 <= z3 ){                                /* 0 <= z1, z2  and z1<=z3*/
              if ( z2 <= z3 ){                               /* 0 <= z1, z2 <= z3*/
                /* everything ok cannot arrive!*/
                ;
              }else{                                         /* 0 <= z1 <= z3 < = z2*/
                /* permuting U2 U3*/
                sign = -sign;
                swap=x2; x2=x3; x3=swap;
                swap=y2; y2=y3; y3=swap;
                swap=z2; z2=z3; z3=swap;
              }}else{                                         /* 0 <= z1, z2  and z3<z1*/
                if ( z1 <= z2 ){                              /* 0 <= z3 < z1 <= z2*/
                  /* permuting U2 U3*/
                  sign = -sign;
                  swap=x2; x2=x3; x3=swap;
                  swap=y2; y2=y3; y3=swap;
                  swap=z2; z2=z3; z3=swap;
                }else{                                        /* 0 <= z2, z3 < z1*/
                  /* permuting U1 U3*/
                  sign = -sign;
                  swap=x1; x1=x3; x3=swap;
                  swap=y1; y1=y3; y3=swap;
                  swap=z1; z1=z3; z3=swap;
                }
              }
            if ((m3=ABDPY_det2x2(x1,y1,x2,y2))==0)
              {                 /* the lattice is twice undifined, z3 can be replaced by 0.0*/
                z3=0.0;
                /* now at least two z values are null*/
                if (sign == 1)
                  if ( z1 == 0.0 )									/* z1 = z3 = 0 <= z2*/
                    if ( z2 == 0.0)
                      return 0;
                    else
                      return   ABDPY_det2x2(x3,y3,x1,y1) ;
                  else								              	/* z2 = z3 = 0 < z1*/
                    return   ABDPY_det2x2(x2,y2,x3,y3) ;
                else
                  if ( z1 == 0.0 )									/* z1 = z3 = 0 <= z2*/
                    if ( z2 == 0.0)
                      return 0;
                    else
                      return   ABDPY_det2x2(x1,y1,x3,y3) ;
                  else								              	/* z2 = z3 = 0 < z1*/
                    return   ABDPY_det2x2(x3,y3,x2,y2) ;
              }
          }
      }
  }







  int ABDPY__det3x3(double& x1, double& y1, double& z1, 
                    double& x2, double& y2, double& z2,
                    double& x3, double& y3, double& z3,
                    short& bb1, short& bb2, short& bb3)
    /* ui refer to vector (xi,yi,zi)*/
    /* bbi is number of bits used in ui + 1023 (IEEE norm)*/
  {
    register short bb;                             /* for lazy evaluation*/
    static   double det;
    static union {double number; short adr[sizeof(double)/sizeof(short)]; } eps;     /* for lazy evaluation*/
    eps.number = 0.0;
    register double swap;
    static   double xx,yy,zz,xxx,yyy/*,zzz*/;
    static   int sign;            /* difference between the sign of the original*/
    /* determinant and the current one (value 1 or -1)*/
    static   int m1,m2,m3;        /* sign of minors*/
    static   double k;            /* quotient in some special cases (null minors)*/
    static   double xv1,yv1,zv1;  /* v1 vector*/
    static   double xv2,yv2,zv2;  /* v2 vector*/
    static   int u_to_v;    	  /* difference between ui and vi*/
    static   double cccx[106];    /* stack of C(l) points*/
    static   double cccy[106];    /* stack of C(l) points*/
    static   double cccz[106];    /* stack of C(l) points*/
    register double* cx;		  /* current place in the stack*/
    register double* cy;		  /* current place in the stack*/
    register double* cz;		  /* current place in the stack*/
    static   double ccx,ccy,ccz;  /* C(l) in step 1.4*/

    sign=1;
    /* first iteration*/
    {
      /* making z coordinates positive and permuting the entries*/
      /* so that z3 is the biggest one*/
      {
        if ( 0.0 <= z1 ){                                   /* 0 <= z1*/
          if ( 0.0 <= z2 ){                                  /* 0 <= z1, z2*/
            if ( z1 <= z3 ){                                /* 0 <= z1, z2  and z1<=z3*/
              if ( z2 <= z3 ){                               /* 0 <= z1, z2 <= z3*/
                /* everything ok*/
                ;
              }else{                                         /* 0 <= z1 <= z3 < = z2*/
                /* permuting U2 U3*/
                sign = -sign;
                swap=x2; x2=x3; x3=swap;
                swap=y2; y2=y3; y3=swap;
                swap=z2; z2=z3; z3=swap;
                bb2=bb3;
              }}else{                                         /* 0 <= z1, z2  and z3<z1*/
                if ( 0.0 <= z3 ) {                               /* 0 <= z2, z3  and z3<z1*/
                  if ( z1 <= z2 ){                              /* 0 <= z3 < z1 <= z2*/
                    /* permuting U2 U3*/
                    sign = -sign;
                    swap=x2; x2=x3; x3=swap;
                    swap=y2; y2=y3; y3=swap;
                    swap=z2; z2=z3; z3=swap;
                    bb2=bb3;
                  }else{                                        /* 0 <= z2, z3 < z1*/
                    /* permuting U1 U3*/
                    sign = -sign;
                    swap=x1; x1=x3; x3=swap;
                    swap=y1; y1=y3; y3=swap;
                    swap=z1; z1=z3; z3=swap;
                    bb1=bb3;
                  }}else{                                        /* z3 < 0 <= z1, z2*/
                    if ( z1 <= -z3 ){                             /* 0 <= z1, z2 and z1 <= -z3*/
                      if ( z2 <= -z3 ){                             /* 0 <= z1, z2 <= -z3*/
                        /* changing U3 sign*/
                        sign = -sign;
                        x3 = -x3; y3 = -y3; z3 = -z3;
                      }else{                                        /* 0 <= z1 <= -z3 < z2*/
                        /* changing U3 sign and permuting U2 U3 */
                        swap=x2; x2=-x3; x3=swap;
                        swap=y2; y2=-y3; y3=swap;
                        swap=z2; z2=-z3; z3=swap;
                        bb2=bb3;
                      }}else{                                        /* 0 <= z1, z2 and -z3 < z1*/
                        if ( z2 <= z1 ){                              /* 0 <= -z3, z2 <= z1*/
                          /* changing U3 sign and permuting U1 U3*/
                          swap=x1; x1=-x3; x3=swap;
                          swap=y1; y1=-y3; y3=swap;
                          swap=z1; z1=-z3; z3=swap;
                          bb1=bb3;
                        }else{                                        /* 0 <= -z3 < z1 < z2*/
                          /* changing U3 sign and permuting U2 U3*/
                          swap=x2; x2=-x3; x3=swap;
                          swap=y2; y2=-y3; y3=swap;
                          swap=z2; z2=-z3; z3=swap;
                          bb2=bb3;
                        }}}}}else{                                        /* z2 < 0 <= z1*/
                          if ( z1 <= z3 ){                                 /* z2 < 0 <= z1 <= z3*/
                            if ( -z2 <= z3 ){                               /* 0 <= z1, -z2 <= z3*/
                              /* changing U2 sign*/
                              sign = -sign;
                              x2 = -x2; y2 = -y2; z2 = -z2;
                            }else{                                          /* 0 <= z1 <= z3 < -z2*/
                              /* changing U2 sign and permuting U2 U3 */
                              swap=-x2; x2=x3; x3=swap;
                              swap=-y2; y2=y3; y3=swap;
                              swap=-z2; z2=z3; z3=swap;
                              bb2=bb3;
                            }}else{                                           /* z2 < 0 <= z1 and z3<z1*/
                              if ( 0.0 <= z3 ) {                                 /* 0<=z1, -z2, z3 and z3<z1*/
                                if ( -z2 <= z1 ){                                /* 0 <= -z2, z3 <= z1*/
                                  /* changing U2 sign and permuting U1 U3 */
                                  x2 = -x2; y2 = -y2; z2 = -z2;
                                  swap=x1; x1=x3; x3=swap;
                                  swap=y1; y1=y3; y3=swap;
                                  swap=z1; z1=z3; z3=swap;
                                  bb1=bb3;
                                }else{                                           /* 0 <= z3 <= z1 < -z2*/
                                  /* changing U2 sign and permuting U2 U3 */
                                  swap=-x2; x2=x3; x3=swap;
                                  swap=-y2; y2=y3; y3=swap;
                                  swap=-z2; z2=z3; z3=swap;
                                  bb2=bb3;
                                }}else{                                           /* 0<=z1, -z2, -z3*/
                                  if ( z1 <= -z3) {                               /* z2 < 0 <= z1 <= -z3*/
                                    if ( z3 <= z2){                                /* 0 <= z1, -z2 <= -z3*/
                                      /* changing U2 U3 sign*/
                                      x2 = -x2; y2 = -y2; z2 = -z2;
                                      x3 = -x3; y3 = -y3; z3 = -z3;
                                    }else{                                         /* 0 <= z1 <= -z3 < -z2*/
                                      /* changing U2 U3 sign and permuting U2 U3 */
                                      sign = -sign;
                                      swap=-x2; x2=-x3; x3=swap;
                                      swap=-y2; y2=-y3; y3=swap;
                                      swap=-z2; z2=-z3; z3=swap;
                                      bb2=bb3;
                                    }}else {                                        /* z2 < 0 <= -z3 < z1*/
                                      if ( z1 <= -z2 ) {                             /* 0 <= -z3 < z1 <= -z2*/
					/* changing U2 U3 sign and permuting U2 U3 */
					sign = -sign;
					swap=-x2; x2=-x3; x3=swap;
					swap=-y2; y2=-y3; y3=swap;
					swap=-z2; z2=-z3; z3=swap;
					bb2=bb3;
                                      }else{                                         /* 0 <= -z3, -z2 < z1*/
					/* changing U2 U3 sign and permuting U1 U3 */
					sign = -sign;
					x2 = -x2; y2 = -y2; z2 = -z2;
					swap=x1; x1=-x3; x3=swap;
					swap=y1; y1=-y3; y3=swap;
					swap=z1; z1=-z3; z3=swap;
					bb1=bb3;
                                      }}}}}}else{                                          /* z1 < 0*/
                                        if ( 0.0 <= z2 ){                                  /* 0 <= -z1, z2*/
                                          if ( -z1 <= z3 ){                                /* 0 <= -z1, z2  and -z1<=z3*/
                                            if ( z2 <= z3 ){                               /* 0 <= -z1, z2 <= z3*/
                                              /* changing U1 sign*/
                                              sign = -sign;
                                              x1 = -x1; y1 = -y1; z1 = -z1;
                                            }else{                                         /* 0 <= -z1 <= z3 < = z2*/
                                              /* changing U1 sign permuting U2 U3*/
                                              x1 = -x1; y1 = -y1; z1 = -z1;
                                              swap=x2; x2=x3; x3=swap;
                                              swap=y2; y2=y3; y3=swap;
                                              swap=z2; z2=z3; z3=swap;
                                              bb2=bb3;
                                            }}else{                                         /* 0 <= -z1, z2  and z3<-z1*/
                                              if ( 0.0 <= z3 ) {                               /* 0 <= z2, z3  and z3<-z1*/
                                                if ( -z1 <= z2 ){                              /* 0 <= z3 < -z1 <= z2*/
                                                  /* changing U1 sign permuting U2 U3*/
                                                  x1 = -x1; y1 = -y1; z1 = -z1;
                                                  swap=x2; x2=x3; x3=swap;
                                                  swap=y2; y2=y3; y3=swap;
                                                  swap=z2; z2=z3; z3=swap;
                                                  bb2=bb3;
                                                }else{                                        /* 0 <= z2, z3 < -z1*/
                                                  /* changing U1 sign permuting U1 U3*/
                                                  swap=-x1; x1=x3; x3=swap;
                                                  swap=-y1; y1=y3; y3=swap;
                                                  swap=-z1; z1=z3; z3=swap;
                                                  bb1=bb3;
                                                }}else{                                        /* z3 < 0 <= -z1, z2*/
                                                  if ( -z1 <= -z3 ){                             /* 0 <=-z1,z2 and -z1 <= -z3*/
                                                    if ( z2 <= -z3 ){                             /* 0 <= -z1, z2 <= -z3*/
                                                      /* changing U1 U3 sign*/
                                                      x1 = -x1; y1 = -y1; z1 = -z1;
                                                      x3 = -x3; y3 = -y3; z3 = -z3;
                                                    }else{                                        /* 0 <= -z1 <= -z3 < z2*/
                                                      /* changing U1 U3 sign and permuting U2 U3 */
                                                      sign = -sign;
                                                      x1 = -x1; y1 = -y1; z1 = -z1;
                                                      swap=x2; x2=-x3; x3=swap;
                                                      swap=y2; y2=-y3; y3=swap;
                                                      swap=z2; z2=-z3; z3=swap;
                                                      bb2=bb3;
                                                    }}else{                                        /* 0 <=-z1,z2 and -z3 < -z1*/
                                                      if ( z2 <= -z1 ){                              /* 0 <= -z3, z2 <= -z1*/
                                                        /* changing U1 U3 sign and permuting U1 U3*/
                                                        sign = -sign;
                                                        swap=-x1; x1=-x3; x3=swap;
                                                        swap=-y1; y1=-y3; y3=swap;
                                                        swap=-z1; z1=-z3; z3=swap;
                                                        bb1=bb3;
                                                      }else{                                        /* 0 <= -z3 < -z1 < z2*/
                                                        /* changing U1 U3 sign and permuting U2 U3*/
                                                        sign = -sign;
                                                        x1 = -x1; y1 = -y1; z1 = -z1;
                                                        swap=x2; x2=-x3; x3=swap;
                                                        swap=y2; y2=-y3; y3=swap;
                                                        swap=z2; z2=-z3; z3=swap;
                                                        bb2=bb3;
                                                      }}}}}else{                                        /* z2 < 0 <= -z1*/
                                                        if ( -z1 <= z3 ){                                 /* z2 < 0 <= -z1 <= z3*/
                                                          if ( -z2 <= z3 ){                               /* 0 <= -z1, -z2 <= z3*/
                                                            /* changing U1U2 sign*/
                                                            x1 = -x1; y1 = -y1; z1 = -z1;
                                                            x2 = -x2; y2 = -y2; z2 = -z2;
                                                          }else{                                          /* 0 <= -z1 <= z3 < -z2*/
                                                            /* changing U1 U2 sign and permuting U2 U3 */
                                                            sign = -sign;
                                                            x1 = -x1; y1 = -y1; z1 = -z1;
                                                            swap=-x2; x2=x3; x3=swap;
                                                            swap=-y2; y2=y3; y3=swap;
                                                            swap=-z2; z2=z3; z3=swap;
                                                            bb2=bb3;
                                                          }}else{                                           /* z2<0 <= -z1 and z3<-z1*/
                                                            if ( 0.0 <= z3 ) {                                 /* 0<=-z1,-z2,z3 and z3<-z1*/
                                                              if ( z2 >= z1 ){                                 /* 0 <= -z2, z3 <= -z1*/
                                                                /* changing U1 U2 sign and permuting U1 U3 */
                                                                sign = -sign;
                                                                x2 = -x2; y2 = -y2; z2 = -z2;
                                                                swap=-x1; x1=x3; x3=swap;
                                                                swap=-y1; y1=y3; y3=swap;
                                                                swap=-z1; z1=z3; z3=swap;
                                                                bb1=bb3;
                                                              }else{                                           /* 0 <= z3 <= -z1 < -z2*/
                                                                /* changing U1 U2 sign and permuting U2 U3 */
                                                                sign = -sign;
                                                                x1 = -x1; y1 = -y1; z1 = -z1;
                                                                swap=-x2; x2=x3; x3=swap;
                                                                swap=-y2; y2=y3; y3=swap;
                                                                swap=-z2; z2=z3; z3=swap;
                                                                bb2=bb3;
                                                              }}else{                                           /* 0<=-z1, -z2, -z3*/
                                                                if ( -z1 <= -z3) {                               /* z2 < 0 <= -z1 <= -z3*/
                                                                  if ( z3 <= z2){                                /* 0 <= -z1, -z2 <= -z3*/
                                                                    /* changing U1 U2 U3 sign*/
                                                                    sign = -sign;
                                                                    x1 = -x1; y1 = -y1; z1 = -z1;
                                                                    x2 = -x2; y2 = -y2; z2 = -z2;
                                                                    x3 = -x3; y3 = -y3; z3 = -z3;
                                                                  }else{                                         /* 0 <= -z1 <= -z3 < -z2*/
                                                                    /* changing U1 U2 U3 sign and permuting U2 U3 */
                                                                    x1 = -x1; y1 = -y1; z1 = -z1;
                                                                    swap=-x2; x2=-x3; x3=swap;
                                                                    swap=-y2; y2=-y3; y3=swap;
                                                                    swap=-z2; z2=-z3; z3=swap;
                                                                    bb2=bb3;
                                                                  }}else {                                        /* z2 < 0 <= -z3 < -z1*/
                                                                    if ( -z1 <= -z2 ) {                             /* 0 <= -z3 < -z1 <= -z2*/
                                                                      /* changing U1 U2 U3 sign and permuting U2 U3 */
                                                                      x1 = -x1; y1 = -y1; z1 = -z1;
                                                                      swap=-x2; x2=-x3; x3=swap;
                                                                      swap=-y2; y2=-y3; y3=swap;
                                                                      swap=-z2; z2=-z3; z3=swap;
                                                                      bb2=bb3;
                                                                    }else{                                         /* 0 <= -z3, -z2 < -z1*/
                                                                      /* changing U1 U2 U3 sign and permuting U1 U3 */
                                                                      x2 = -x2; y2 = -y2; z2 = -z2;
                                                                      swap=-x1; x1=-x3; x3=swap;
                                                                      swap=-y1; y1=-y3; y3=swap;
                                                                      swap=-z1; z1=-z3; z3=swap;
                                                                      bb1=bb3;
                                                                    }}}}}}}
      if ((m3=ABDPY_det2x2(x1,y1,x2,y2))==0)
        {                  /* the lattice is undefined, but z3 can be replaced by 0.0*/
          z3=0.0;
          /* redo the permutation*/
          if ( z1 <= z3 ){                                /* 0 <= z1, z2  and z1<=z3*/
            if ( z2 <= z3 ){                               /* 0 <= z1, z2 <= z3*/
              /* everything ok cannot arrive!*/
              ;
            }else{                                         /* 0 <= z1 <= z3 < = z2*/
              /* permuting U2 U3*/
              sign = -sign;
              swap=x2; x2=x3; x3=swap;
              swap=y2; y2=y3; y3=swap;
              swap=z2; z2=z3; z3=swap;
            }}else{                                         /* 0 <= z1, z2  and z3<z1*/
              if ( z1 <= z2 ){                              /* 0 <= z3 < z1 <= z2*/
                /* permuting U2 U3*/
                sign = -sign;
                swap=x2; x2=x3; x3=swap;
                swap=y2; y2=y3; y3=swap;
                swap=z2; z2=z3; z3=swap;
              }else{                                        /* 0 <= z2, z3 < z1*/
                /* permuting U1 U3*/
                sign = -sign;
                swap=x1; x1=x3; x3=swap;
                swap=y1; y1=y3; y3=swap;
                swap=z1; z1=z3; z3=swap;
              }
            }
          if ((m3=ABDPY_det2x2(x1,y1,x2,y2))==0)
            {                 /* the lattice is twice undifined, z3 can be replaced by 0.0*/
              z3=0.0;
              /* now at least two z values are null*/
              if (sign == 1)
                if ( z1 == 0.0 )									/* z1 = z3 = 0 <= z2*/
                  if ( z2 == 0.0)
                    return 0;
                  else
                    return   ABDPY_det2x2(x3,y3,x1,y1) ;
                else								              	/* z2 = z3 = 0 < z1*/
                  return   ABDPY_det2x2(x2,y2,x3,y3) ;
              else
                if ( z1 == 0.0 )									/* z1 = z3 = 0 <= z2*/
                  if ( z2 == 0.0)
                    return 0;
                  else
                    return   ABDPY_det2x2(x1,y1,x3,y3) ;
                else								              	/* z2 = z3 = 0 < z1*/
                  return   ABDPY_det2x2(x3,y3,x2,y2) ;
            }
          {/*LAZY EVALUATION*/
            bb1 = ((ABDPY_shortdouble(x1))& 0x7ff0)>>4;
            bb = ((ABDPY_shortdouble(y1))& 0x7ff0)>>4; if (bb>bb1) bb1=bb;
            bb = ((ABDPY_shortdouble(z1))& 0x7ff0)>>4; if (bb>bb1) bb1=bb;
            bb2 = ((ABDPY_shortdouble(x2))& 0x7ff0)>>4;
            bb = ((ABDPY_shortdouble(y2))& 0x7ff0)>>4; if (bb>bb2) bb2=bb;
            bb = ((ABDPY_shortdouble(z2))& 0x7ff0)>>4; if (bb>bb2) bb2=bb;
            bb3 = ((ABDPY_shortdouble(x3))& 0x7ff0)>>4;
            bb = ((ABDPY_shortdouble(y3))& 0x7ff0)>>4; if (bb>bb3) bb3=bb;
            bb = ((ABDPY_shortdouble(z3))& 0x7ff0)>>4; if (bb>bb3) bb3=bb;
            bb = bb1+bb2+bb3 -2094;
            if (bb < 64) return 0;
            det = x1*((y2*z3)-(y3*z2)) - y1*((x2*z3)-(x3*z2)) + z1*((x2*y3)-(x3*y2));
            *eps.adr = (bb<<4); 
            if (det >  eps.number) return  sign;
            if (det < -eps.number) return -sign;
            if (bb<1023) return 0;
          }
        }
    }	
    while(1)
      {
        if (z3 == 0.0) return 0;                            /* z1=z2=z3=0.0*/
        /* compute the remaining minors*/
        m1 = ABDPY_det2x2(x2,y2,x3,y3);
        m2 = ABDPY_det2x2(x3,y3,x1,y1);
        /* if the three minors have the same sign*/
        if ( ( m1 >= 0 ) && ( m2 >= 0) && (m3 >= 0)) return sign ;
        if ( ( m1 <= 0 ) && ( m2 <= 0) && (m3 <= 0)) return -sign ;
        /* if one minor is null the corresponding z value can be replaced by 0.0*/
        /* m3==0 already tested*/
        if ( m1 == 0 )
          {
            z1 = 0.0;
            if (x2 != 0.0)
              {
		k = floor(x3/x2);    /* as in 2D case k can be computed directly*/
		x3 -=k*x2;
		y3 -=k*y2;
		z3 -=k*z2;
              }
            else       /* x2=0.0*/
              {
		k = floor(y3/y2);    /* as in 2D case k can be computed directly*/
		y3 -=k*y2;
		z3 -=k*z2;
              }
          }
        else if ( m2 == 0 )            /* m1 != 0*/
          {
            z2 = 0.0;
            if (x1 != 0.0)
              {
		k = floor(x3/x1);    /* as in 2D case k can be computed directly*/
		x3 -=k*x1;
		y3 -=k*y1;
		z3 -=k*z1;
              }
            else       /* x1=0.0*/
              {
		k = floor(y3/y1);    /* as in 2D case k can be computed directly*/
		y3 -=k*y1;
		z3 -=k*z1;
              }
          }
        else                           /* m1 != 0 and m2 != 0 and m3 != 0*/
          {													/* step 1 : compute R*/
            /* substep 1.1*/
            /* computing v1 v2 such that v1 + v2 and u3 is in the positive wedge v1 v2*/
            if ( m3 > 0 )                   /*           m3>0*/
              if ( m1 > 0 )                  /* m1>0      m3>0*/
                if ( m2 > 0 )                 /* m1>0 m2>0 m3>0*/
                  {
                    u_to_v = 1;
                    xv1 = -x1; yv1 = -y1; zv1 = -z1;
                    xv2 = -x2; yv2 = -y2; zv2 = -z2;
                  }
                else                          /* m1>0 m2<0 m3>0*/
                  {
                    u_to_v = 2;
                    xv1 =  x2; yv1 =  y2; zv1 =  z2;
                    xv2 = -x1; yv2 = -y1; zv2 = -z1;
                  }
              else                           /* m1<0      m3>0*/
                if ( m2 > 0 )                 /* m1<0 m2>0 m3>0*/
                  {
                    u_to_v = 3;
                    xv1 = -x2; yv1 = -y2; zv1 = -z2;
                    xv2 =  x1; yv2 =  y1; zv2 =  z1;
                  }
                else                          /* m1<0 m2<0 m3>0*/
                  {
                    u_to_v = 4;
                    xv1 =  x1; yv1 =  y1; zv1 =  z1;
                    xv2 =  x2; yv2 =  y2; zv2 =  z2;
                  }
            else                            /*           m3<0*/
              if ( m1 > 0 )                  /* m1>0      m3<0*/
                if ( m2 > 0 )                 /* m1>0 m2>0 m3<0*/
                  {
                    u_to_v = 5;
                    xv1 =  x2; yv1 =  y2; zv1 =  z2;
                    xv2 =  x1; yv2 =  y1; zv2 =  z1;
                  }
                else                          /* m1>0 m2<0 m3<0*/
                  {
                    u_to_v = 6;
                    xv1 =  x1; yv1 =  y1; zv1 =  z1;
                    xv2 = -x2; yv2 = -y2; zv2 = -z2;
                  }
              else                           /* m1<0      m3<0*/
                if ( m2 > 0 )                 /* m1<0 m2>0 m3<0*/
                  {
                    u_to_v = 7;
                    xv1 = -x1; yv1 = -y1; zv1 = -z1;
                    xv2 =  x2; yv2 =  y2; zv2 =  z2;
                  }
                else                          /* m1<0 m2<0 m3<0*/
                  {
                    u_to_v = 8;
                    xv1 = -x2; yv1 = -y2; zv1 = -z2;
                    xv2 = -x1; yv2 = -y1; zv2 = -z1;
                  }
            /* substep 1.2*/
            /* verify if u3 is in the base cell of lattice v1 v2*/
            if ((  ABDPY_det2x2(xv1, yv1, x3-xv2, y3-yv2) > 0 ) ||
                (  ABDPY_det2x2(x3-xv1, y3-yv1, xv2, yv2) > 0 ))
              {
                if (  ABDPY_det2x2(x3, y3, xv1+xv2, yv1+yv2) > 0 )
                  /* u3 in the wedge v1 v1+v2*/
                  {
                    /* substep 1.3*/
                    *(cx=cccx) = xv1;
                    *(cy=cccy) = yv1;
                    *(cz=cccz) = zv1;
                    while (1)
                      { 
                        zz=*cz+zv2;
                        if ((*cz<0)&&(zz<0)) return (m3>0.0) ? sign : -sign ;
                        if ((*cz>ABDPY_2EXP53)&&(zz>ABDPY_2EXP53))
                          return (m3>0.0) ? -sign : sign;
                        xx=*cx+*cx; 
                        yy=*cy+*cy; 
                        if ( ABDPY_det2x2(xv2, yv2, x3-xx, y3-yy) >= 0)
                          break; /* goto 1.4*/
                        xxx = xx + xv2;
                        yyy = yy + yv2;
                        if ( ABDPY_det2x2(x3, y3, xxx, yyy) >= 0)
                          { *++cx =  xx; *++cy =  yy; zz= *cz+*cz; *++cz = zz;}
                        else
                          { *++cx = xxx; *++cy = yyy; zz+= *cz;   *++cz = zz;}
                      }
                    /* substep 1.4*/
                    ccx = *cx;
                    ccy = *cy;
                    ccz = *cz;
                    while (--cx >= cccx)
                      { 
                        --cy;--cz;
                        zz=ccz+zv2;
                        if ((ccz<0)&&(zz<0)) return (m3>0.0) ? sign : -sign ;
                        if ((ccz>ABDPY_2EXP53)&&(zz>ABDPY_2EXP53))
                          return (m3>0.0) ? -sign : sign;
                        xx=ccx+*cx; 
                        yy=ccy+*cy; 
                        if ( ABDPY_det2x2(xv2, yv2, x3-xx, y3-yy) >= 0)
                          continue;
                        xxx = xx + xv2;
                        yyy = yy + yv2;
                        if ( ABDPY_det2x2(x3, y3, xxx, yyy) >= 0)
                          { ccx =  xx; ccy =  yy; ccz+= *cz;}
                        else
                          { ccx = xxx; ccy = yyy; ccz = zz+*cz;}
                      }
                    /* substep 1.5*/
                    xxx = ccx + xv2;
                    yyy = ccy + yv2;
                    if ( ABDPY_det2x2(xv1, yv1, x3-xxx, y3-yyy) > 0)
                      { ccx = xxx; ccy = yyy; ccz += zv2; }
                  }
                else
                  /* u3 in the wedge v1+v2 v2*/
                  {
                    /* substep 1.3*/
                    *(cx=cccx) = xv2;
                    *(cy=cccy) = yv2;
                    *(cz=cccz) = zv2;
                    while (1)
                      { 
                        zz=*cz+zv1;
                        if ((*cz<0)&&(zz<0)) return (m3>0.0) ? sign : -sign ;
                        if ((*cz>ABDPY_2EXP53)&&(zz>ABDPY_2EXP53))
                          return (m3>0.0) ? -sign : sign;
                        xx=*cx+*cx; 
                        yy=*cy+*cy; 
                        if ( ABDPY_det2x2(xv1, yv1, x3-xx, y3-yy) <= 0)
                          break; /* goto 1.4*/
                        xxx = xx + xv1;
                        yyy = yy + yv1;
                        if ( ABDPY_det2x2(x3, y3, xxx, yyy) <= 0)
                          { *++cx =  xx; *++cy =  yy; zz= *cz+*cz; *++cz = zz;}
                        else
                          { *++cx = xxx; *++cy = yyy; zz+= *cz; *++cz = zz;}
                      }
                    /* substep 1.4*/
                    ccx = *cx;
                    ccy = *cy;
                    ccz = *cz;
                    while (--cx >= cccx)
                      { 
                        --cy;--cz;
                        zz=ccz+zv1;
                        if ((ccz<0)&&(zz<0)) return (m3>0.0) ? sign : -sign ;
                        if ((ccz>ABDPY_2EXP53)&&(zz>ABDPY_2EXP53))
                          return (m3>0.0) ? -sign : sign;
                        xx=ccx+*cx; 
                        yy=ccy+*cy; 
                        if ( ABDPY_det2x2(xv1, yv1, x3-xx, y3-yy) <= 0)
                          continue;
                        xxx = xx + xv1;
                        yyy = yy + yv1;
                        if ( ABDPY_det2x2(x3, y3, xxx, yyy) <= 0)
                          { ccx =  xx; ccy =  yy; ccz+= *cz;}
                        else
                          { ccx = xxx; ccy = yyy; ccz = zz+*cz;}
                      }
                    /* substep 1.5*/
                    xxx = ccx + xv1;
                    yyy = ccy + yv1;
                    if ( ABDPY_det2x2(xv2, yv2, x3-xxx, y3-yyy) < 0)
                      { ccx = xxx; ccy = yyy; ccz += zv1; }
                  }	  
                switch (u_to_v)
                  {
                  case 1:
                  case 8:
                    ccx += xv1+xv2; ccy += yv1+yv2; ccz += zv1+zv2; break;
                  case 3:
                  case 7:
                    ccx += xv1; ccy += yv1; ccz += zv1; break;
                  case 2:
                  case 6:
                    ccx += xv2; ccy += yv2; ccz += zv2; break;
                  }
                x3 -= ccx; y3 -= ccy; z3 -= ccz;
              }
            else             /* u3 is originally in the base cell of lattice v1 v2*/
              switch (u_to_v)
                {
                case 1:
                case 8:
                  x3 -= xv1+xv2; y3 -= yv1+yv2; z3 -= zv1+zv2; break;
                case 3:
                case 7:
                  x3 -= xv1; y3 -= yv1; z3 -= zv1; break;
                case 2:
                case 6:
                  x3 -= xv2; y3 -= yv2; z3 -= zv2; break;
                }
          }
        /* substep 1.6*/
        if (z3<0.0) return (m3>0.0) ? -sign : sign;
        if (z3>z1+z2) return (m3>0.0) ? sign : -sign ;
        {
          /* step 2*/
          /* finding R' with z >= 0*/
          if ( ABDPY_det2x2(x2-x1, y2-y1, x3-x1, y3-y1) > 0 )
            {if ( ABDPY_det2x2(x3, y3, x1+x2, y1+y2) > 0 )
              {if ( m3 > 0.0 )                                              /* 0 or u1*/
                {if ( z1 > z2 )
                  if ( z3 > z1 ) return sign;
                  else {if ( z3 + z3 > z1 ) {sign=-sign;z3=z1-z3; x3=x1-x3; y3=y1-y3;}}
                else
                  if ( (ccz=z3 + z3) > z1 + z2 ) return sign;
                  else {if ( z3 > z1 )  {z3-=z1; x3-=x1; y3-=y1;}
                  else if ( ccz > z1 ) {sign=-sign;z3=z1-z3; x3=x1-x3; y3=y1-y3;}}
                }
              else                                                         /* u1+u2 or u2*/
                {if ( z2 < z1 )
                  if ( z3 < z2 ) return sign;
                  else {z3-=z2;x3-=x2;y3-=y2;
                  if ( z3+z3 > z1 ) {sign=-sign;z3=z1-z3; x3=x1-x3; y3=y1-y3;}}
                else
                  if ( z3 + z3 < z1 + z2 ) return sign;
                  else {if ( z3 < z2 ) {sign=-sign;z3=z2-z3; x3=x2-x3; y3=y2-y3;}
                  else{z3-=z2;x3-=x2;y3-=y2;
                  if ( z3+z3 > z1 ) {sign=-sign;z3=z1-z3; x3=x1-x3; y3=y1-y3;}}}
                }
              }
            else
              {if ( m3 > 0.0 )                                              /* 0 or u2*/
                {if ( z2 > z1 )
                  if ( z3 > z2 ) return sign;
                  else {if ( z3 + z3 > z2 ) {sign=-sign;z3=z2-z3; x3=x2-x3; y3=y2-y3;}}
                else
                  if ( (ccz=z3 + z3) > z1 + z2 ) return sign;
                  else {if ( z3 > z2 ) {z3-=z2;x3-=x2;y3-=y2;}
                  else if ( ccz > z2 ) {sign=-sign;z3=z2-z3; x3=x2-x3; y3=y2-y3;}}
                }
              else                                                         /* u1+u2 or u1*/
                {if ( z1 < z2 )
                  if ( z3 < z1 ) return sign;
                  else {z3-=z1;x3-=x1;y3-=y1;
                  if ( z3+z3 > z2 ) {sign=-sign;z3=z2-z3; x3=x2-x3; y3=y2-y3;}}
                else
                  if ( z3 + z3 < z1 + z2 ) return sign;
                  else {if ( z3 < z1 ) {sign=-sign;z3=z1-z3; x3=x1-x3; y3=y1-y3;}
                  else{z3-=z1;x3-=x1;y3-=y1;
                  if ( z3+z3 > z2 ) {sign=-sign;z3=z2-z3; x3=x2-x3; y3=y2-y3;}}}
                }
              }
            }
          else
            {if ( ABDPY_det2x2(x3, y3, x1+x2, y1+y2) > 0 )
              {if ( m3 > 0.0 )                                              /* u1+u2 or u1*/
                {if ( z1 < z2 )
                  if ( z3 < z1 ) return -sign;
                  else {z3-=z1;x3-=x1;y3-=y1;
                  if ( z3+z3 > z2 ) {sign=-sign;z3=z2-z3; x3=x2-x3; y3=y2-y3;}}
                else
                  if ( z3 + z3 < z1 + z2 ) return -sign;
                  else {if ( z3 < z1 ) {sign=-sign;z3=z1-z3; x3=x1-x3; y3=y1-y3;}
                  else{z3-=z1;x3-=x1;y3-=y1;
                  if ( z3+z3 > z2 ) {sign=-sign;z3=z2-z3; x3=x2-x3; y3=y2-y3;}}}
                }
              else                                                         /* 0 or u2*/
                {if ( z2 > z1 )
                  if ( z3 > z2 ) return -sign;
                  else {if ( z3 + z3 > z2 ) {sign=-sign;z3=z2-z3; x3=x2-x3; y3=y2-y3;}}
                else
                  if ( (ccz=z3 + z3) > z1 + z2 ) return -sign;
                  else {if ( z3 > z2 ) {z3-=z2;x3-=x2;y3-=y2;}
                  else if ( ccz > z2 ) {sign=-sign;z3=z2-z3; x3=x2-x3; y3=y2-y3;}}
                }
              }
            else
              {if ( m3 > 0.0 )                                              /* u1+u2 or u2*/
                {if ( z2 < z1 )
                  if ( z3 < z2 ) return -sign;
                  else {z3-=z2;x3-=x2;y3-=y2;
                  if ( z3+z3 > z1 ) {sign=-sign;z3=z1-z3; x3=x1-x3; y3=y1-y3;}}
                else
                  if ( z3 + z3 < z1 + z2 ) return -sign;
                  else {if ( z3 < z2 ) {sign=-sign;z3=z2-z3; x3=x2-x3; y3=y2-y3;}
                  else{z3-=z2;x3-=x2;y3-=y2;
                  if ( z3+z3 > z1 ) {sign=-sign;z3=z1-z3; x3=x1-x3; y3=y1-y3;}}}
                }
              else                                                         /* 0 or u1*/
                {if ( z1 > z2 )
                  if ( z3 > z1 ) return -sign;
                  else {if ( z3 + z3 > z1 ) {sign=-sign;z3=z1-z3; x3=x1-x3; y3=y1-y3;}}
                else
                  if ( (ccz=z3 + z3) > z1 + z2 ) return -sign;
                  else {if ( z3 > z1 )  {z3-=z1; x3-=x1; y3-=y1;}
                  else if ( ccz > z1 ) {sign=-sign;z3=z1-z3; x3=x1-x3; y3=y1-y3;}}
                }
              }
            }
        } /* end R'*/
        {/*LAZY EVALUATION*/
          bb3 = ((ABDPY_shortdouble(x3))& 0x7ff0)>>4;
          bb = ((ABDPY_shortdouble(y3))& 0x7ff0)>>4; if (bb>bb3) bb3=bb;
          bb = ((ABDPY_shortdouble(z3))& 0x7ff0)>>4; if (bb>bb3) bb3=bb;
          bb = bb1+bb2+bb3 -2094;
          if (bb < 64) return 0;
          det = x1*((y2*z3)-(y3*z2)) - y1*((x2*z3)-(x3*z2)) + z1*((x2*y3)-(x3*y2));
          *eps.adr = (bb<<4);
          if (det >  eps.number) return  sign;
          if (det < -eps.number) return -sign;
          if (bb<1023) return 0;
        }
        /* further iteration  0 <= z1,z2,z3 <= max(z1,z2)/2*/
        {
          /* making z coordinates positive and permuting the entries*/
          /* so that z3 is the biggest one*/
          if ( z1 <= z3 ){                                /* 0 <= z1, z2  and z1<=z3*/
            if ( z2 <= z3 ){                               /* 0 <= z1, z2 <= z3*/
              /* everything ok cannot arrive!*/
              ;
            }else{                                         /* 0 <= z1 <= z3 < = z2*/
              /* permuting U2 U3*/
              sign = -sign;
              swap=x2; x2=x3; x3=swap;
              swap=y2; y2=y3; y3=swap;
              swap=z2; z2=z3; z3=swap;
              bb2=bb3;
            }}else{                                         /* 0 <= z1, z2  and z3<z1*/
              if ( z1 <= z2 ){                              /* 0 <= z3 < z1 <= z2*/
                /* permuting U2 U3*/
                sign = -sign;
                swap=x2; x2=x3; x3=swap;
                swap=y2; y2=y3; y3=swap;
                swap=z2; z2=z3; z3=swap;
                bb2=bb3;
              }else{                                        /* 0 <= z2, z3 < z1*/
                /* permuting U1 U3*/
                sign = -sign;
                swap=x1; x1=x3; x3=swap;
                swap=y1; y1=y3; y3=swap;
                swap=z1; z1=z3; z3=swap;
                bb1=bb3;
              }
            }}
        if ((m3=ABDPY_det2x2(x1,y1,x2,y2))==0)
          {                  /* the lattice is undefined, but z3 can be replaced by 0.0*/
            z3=0.0;
            /* redo the permutation*/
            if ( z1 <= z3 ){                                /* 0 <= z1, z2  and z1<=z3*/
              if ( z2 <= z3 ){                               /* 0 <= z1, z2 <= z3*/
                /* everything ok cannot arrive!*/
                ;
              }else{                                         /* 0 <= z1 <= z3 < = z2*/
                /* permuting U2 U3*/
                sign = -sign;
                swap=x2; x2=x3; x3=swap;
                swap=y2; y2=y3; y3=swap;
                swap=z2; z2=z3; z3=swap;
              }}else{                                         /* 0 <= z1, z2  and z3<z1*/
                if ( z1 <= z2 ){                              /* 0 <= z3 < z1 <= z2*/
                  /* permuting U2 U3*/
                  sign = -sign;
                  swap=x2; x2=x3; x3=swap;
                  swap=y2; y2=y3; y3=swap;
                  swap=z2; z2=z3; z3=swap;
                }else{                                        /* 0 <= z2, z3 < z1*/
                  /* permuting U1 U3*/
                  sign = -sign;
                  swap=x1; x1=x3; x3=swap;
                  swap=y1; y1=y3; y3=swap;
                  swap=z1; z1=z3; z3=swap;
                }
              }
            if ((m3=ABDPY_det2x2(x1,y1,x2,y2))==0)
              {                 /* the lattice is twice undifined, z3 can be replaced by 0.0*/
                z3=0.0;
                /* now at least two z values are null*/
                if (sign == 1)
                  if ( z1 == 0.0 )									/* z1 = z3 = 0 <= z2*/
                    if ( z2 == 0.0)
                      return 0;
                    else
                      return   ABDPY_det2x2(x3,y3,x1,y1) ;
                  else								              	/* z2 = z3 = 0 < z1*/
                    return   ABDPY_det2x2(x2,y2,x3,y3) ;
                else
                  if ( z1 == 0.0 )									/* z1 = z3 = 0 <= z2*/
                    if ( z2 == 0.0)
                      return 0;
                    else
                      return   ABDPY_det2x2(x1,y1,x3,y3) ;
                  else								              	/* z2 = z3 = 0 < z1*/
                    return   ABDPY_det2x2(x3,y3,x2,y2) ;
              }
            {/*LAZY EVALUATION*/
              bb1 = ((ABDPY_shortdouble(x1))& 0x7ff0)>>4;
              bb = ((ABDPY_shortdouble(y1))& 0x7ff0)>>4; if (bb>bb1) bb1=bb;
              bb = ((ABDPY_shortdouble(z1))& 0x7ff0)>>4; if (bb>bb1) bb1=bb;
              bb2 = ((ABDPY_shortdouble(x2))& 0x7ff0)>>4;
              bb = ((ABDPY_shortdouble(y2))& 0x7ff0)>>4; if (bb>bb2) bb2=bb;
              bb = ((ABDPY_shortdouble(z2))& 0x7ff0)>>4; if (bb>bb2) bb2=bb;
              bb3 = ((ABDPY_shortdouble(x3))& 0x7ff0)>>4;
              bb = ((ABDPY_shortdouble(y3))& 0x7ff0)>>4; if (bb>bb3) bb3=bb;
              bb = ((ABDPY_shortdouble(z3))& 0x7ff0)>>4; if (bb>bb3) bb3=bb;
              bb = bb1+bb2+bb3 -2094;
              if (bb < 64) return 0;
              det = x1*((y2*z3)-(y3*z2)) - y1*((x2*z3)-(x3*z2)) + z1*((x2*y3)-(x3*y2));
              *eps.adr = (bb<<4);
              if (det >  eps.number) return  sign;
              if (det < -eps.number) return -sign;
              if (bb<1023) return 0;
            }
          }
      }
  }

  int ABDPY_det3x3(double x1, double y1, double z1, 
                   double x2, double y2, double z2,
                   double x3, double y3, double z3)
  {
    short bb1,bb2,bb3;                           /* for lazy evaluation*/
    {
      register short bb;                           /* for lazy evaluation*/
      static double det;
      static union {double number; short adr[sizeof(double)/sizeof(short)]; } eps;     /* for lazy evaluation*/
      eps.number = 0.0;
      static int sign;
      sign=1;

      {/*LAZY EVALUATION*/
        bb1 = ((ABDPY_shortdouble(x1))& 0x7ff0)>>4;
        bb = ((ABDPY_shortdouble(y1))& 0x7ff0)>>4; if (bb>bb1) bb1=bb;
        bb = ((ABDPY_shortdouble(z1))& 0x7ff0)>>4; if (bb>bb1) bb1=bb;
        bb2 = ((ABDPY_shortdouble(x2))& 0x7ff0)>>4;
        bb = ((ABDPY_shortdouble(y2))& 0x7ff0)>>4; if (bb>bb2) bb2=bb;
        bb = ((ABDPY_shortdouble(z2))& 0x7ff0)>>4; if (bb>bb2) bb2=bb;
        bb3 = ((ABDPY_shortdouble(x3))& 0x7ff0)>>4;
        bb = ((ABDPY_shortdouble(y3))& 0x7ff0)>>4; if (bb>bb3) bb3=bb;
        bb = ((ABDPY_shortdouble(z3))& 0x7ff0)>>4; if (bb>bb3) bb3=bb;
        bb = bb1+bb2+bb3 -2094;
        det = x1*((y2*z3)-(y3*z2)) - y1*((x2*z3)-(x3*z2)) + z1*((x2*y3)-(x3*y2));
	*eps.adr = (bb<<4); 
	if (det >  eps.number) return  sign;
	if (det < -eps.number) return -sign;
        if (bb<1023) return 0;
      }
    }
    return ABDPY__det3x3(x1,y1,z1,x2,y2,z2,x3,y3,z3,bb1,bb2,bb3);
  }



  int ABDPY_secure_det2x2(double x1, double y1, double x2, double y2)
  {
    double p,q;
    short bb1,bb;
    {/* TEST IF DATA ARE CORRECT*/
      bb1 = ((ABDPY_shortdouble(x1))& 0x7ff0)>>4;
      bb = ((ABDPY_shortdouble(y1))& 0x7ff0)>>4; if (bb>bb1) bb1=bb;
      bb = ((ABDPY_shortdouble(x2))& 0x7ff0)>>4; if (bb>bb1) bb1=bb;
      bb = ((ABDPY_shortdouble(y2))& 0x7ff0)>>4; if (bb>bb1) bb1=bb;
      if (bb1>53+1023)  { return 2;}
      if (x1!=floor(x1)){ return 3;}
      if (y1!=floor(y1)){ return 3;}
      if (x2!=floor(x2)){ return 3;}
      if (y2!=floor(y2)){ return 3;}
    }
    /* LAZY EVALUATION*/
    p =x1*y2;
    q =x2*y1;
    if (p > q) return 1;
    if (p < q) return -1;
    return ABDPY__det2x2(x1,y1,x2,y2);
  }


  int ABDPY_secure_det3x3(double x1, double y1, double z1, 
                          double x2, double y2, double z2,
                          double x3, double y3, double z3)
  {
    short bb1,bb2,bb3;                           /* for lazy evaluation*/
    {
      register short bb;                           /* for lazy evaluation*/
      static double det;
      static union {double number; short adr[sizeof(double)/sizeof(short)]; } eps;     /* for lazy evaluation*/
      eps.number = 0.0;
      static int sign;
      sign=1;

      {/*LAZY EVALUATION*/
        /* log2(A) = (((ABDPY_shortdouble(A)) & 0x7ff0)>>4)-1023; for A!= 0.0*/
        bb1 = ((ABDPY_shortdouble(x1))& 0x7ff0)>>4;
        bb = ((ABDPY_shortdouble(y1))& 0x7ff0)>>4; if (bb>bb1) bb1=bb;
        bb = ((ABDPY_shortdouble(z1))& 0x7ff0)>>4; if (bb>bb1) bb1=bb;
        bb2 = ((ABDPY_shortdouble(x2))& 0x7ff0)>>4;
        bb = ((ABDPY_shortdouble(y2))& 0x7ff0)>>4; if (bb>bb2) bb2=bb;
        bb = ((ABDPY_shortdouble(z2))& 0x7ff0)>>4; if (bb>bb2) bb2=bb;
        bb3 = ((ABDPY_shortdouble(x3))& 0x7ff0)>>4;
        bb = ((ABDPY_shortdouble(y3))& 0x7ff0)>>4; if (bb>bb3) bb3=bb;
        bb = ((ABDPY_shortdouble(z3))& 0x7ff0)>>4; if (bb>bb3) bb3=bb;
        {/* TEST IF DATA ARE CORRECT*/
          if (bb1>51+1023)  { return 2;}
          if (bb2>51+1023)  { return 2;}
          if (bb3>51+1023)  { return 2;}
          if (x1!=floor(x1)){ return 3;}
          if (y1!=floor(y1)){ return 3;}
          if (z1!=floor(z1)){ return 3;}
          if (x2!=floor(x2)){ return 3;}
          if (y2!=floor(y2)){ return 3;}
          if (z2!=floor(z2)){ return 3;}
          if (x3!=floor(x3)){ return 3;}
          if (y3!=floor(y3)){ return 3;}
          if (z3!=floor(z3)){ return 3;}
        }
        bb = bb1+bb2+bb3 -2094;
        /* if bi is the number of bits used to represent Ui coordinates*/
        /* and if the arithmetic used ba bits for the mantissa*/
        /* then the the determinant computed*/
        /* by the floating point arithmetic is exact with a precision*/
        /* of 2^b     where b=b1+b2+b3-ba+5*/
        /* with IEEE arithmetic, the exponent stored is bb = 1023+b*/
        /* thus bb = 1023 + (bb1+bb2+bb3-3*1023) - 53 + 5=bb1+bb2+bb3 -2094;*/
        det = x1*((y2*z3)-(y3*z2)) - y1*((x2*z3)-(x3*z2)) + z1*((x2*y3)-(x3*y2));
        /* compute epsilon = 2^b  using IEEE format for double*/
	*eps.adr = (bb<<4); 
	if (det >  eps.number) return  sign;
	if (det < -eps.number) return -sign;
        if (bb<1023) return 0;    /* the computaion have been done in full precision*/
      }
    }
    return ABDPY__det3x3(x1,y1,z1,x2,y2,z2,x3,y3,z3,bb1,bb2,bb3);
  }

  // ---------- Wrap

  int idet2x2_sign(double a, double b, double c, double d) {
    return ABDPY_det2x2(a,b,
                        c,d);
  }

  int not_lazy_idet2x2_sign(double a, double b, double c, double d) {
    return ABDPY__not_lazy_det2x2(a,b,
                                  c,d);
  }

  int checked_idet2x2_sign(double a, double b, double c, double d) {
    return ABDPY_secure_det2x2(a,b,
                               c,d);
  }
  
  int idet3x3_sign(double a, double b, double c, double d, double e, double f, double g, double h, double i) {
    return ABDPY_det3x3(a,b,c,
                        d,e,f,
                        g,h,i);
  }

  int not_lazy_idet3x3_sign(double a, double b, double c, double d, double e, double f, double g, double h, double i) {
    return ABDPY__not_lazy_det3x3(a,b,c,
                                  d,e,f,
                                  g,h,i);
  }

  int checked_idet3x3_sign(double a, double b, double c, double d, double e, double f, double g, double h, double i) {
    return ABDPY_secure_det3x3(a,b,c,
                               d,e,f,
                               g,h,i);
  }

}
