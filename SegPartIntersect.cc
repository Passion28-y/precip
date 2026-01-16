/****************************************************************************
 * 
 *      Module:       SegPartIntersect.c
 *
 *      Description:  Contains routines used for maintaining and
 *                    searching the list of inclusion-intersecting
 *                    segments created during the force calculations
 *                    for use with mobility functions that differentiate
 *                    between mobilities inside and outside inclusions.
 *
 *      Includes public functions:
 *          GetIndexFromInclusionID()
 *          FindNodePartIntersects()
 *          SegPartListAppend()
 *          SegPartListClear()
 *          SegPartListInsert
 *          SegPartListLookup()
 *          SegPartListSort()
 *          SegPartListUnsortedLookup()
 *          SegPartListUpdate()
 *
 *      Includes local fucntions()
 *          FindSegPartIntersects()
 *          SegPartListCmp()
 *
 ****************************************************************************/

#ifdef ESHELBY

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Util.h"
#include "Mobility.h"
#include "M33.h"
#include "M44.h"
#include "Topology.h"
#include "Node.h"

int SnapPointToEllipsoidOnPlane(Home_t *home, EInclusion_t *inc,
                                       const real8 planePt[3], const real8 planeN_in[3],
                                       const real8 pin[3], real8 pout[3])
{
    Param_t *param = home->param;

    // Normalize plane normal
    real8 n[3] = { planeN_in[X], planeN_in[Y], planeN_in[Z] };
    real8 nMag = sqrt(n[X]*n[X] + n[Y]*n[Y] + n[Z]*n[Z]);
    if (nMag < 1.0e-20) return 0;
    n[X] /= nMag; n[Y] /= nMag; n[Z] /= nMag;

    // Build orthonormal basis (e1,e2) on the plane
    real8 e1[3];
    if (fabs(n[X]) < 0.9) { e1[X]=1.0; e1[Y]=0.0; e1[Z]=0.0; }
    else                  { e1[X]=0.0; e1[Y]=1.0; e1[Z]=0.0; }
    // e1 = e1 - (e1·n) n
    real8 t = e1[X]*n[X] + e1[Y]*n[Y] + e1[Z]*n[Z];
    e1[X] -= t*n[X]; e1[Y] -= t*n[Y]; e1[Z] -= t*n[Z];
    real8 e1Mag = sqrt(e1[X]*e1[X] + e1[Y]*e1[Y] + e1[Z]*e1[Z]);
    if (e1Mag < 1.0e-20) return 0;
    e1[X] /= e1Mag; e1[Y] /= e1Mag; e1[Z] /= e1Mag;

    real8 e2[3];
    e2[X] = n[Y]*e1[Z] - n[Z]*e1[Y];
    e2[Y] = n[Z]*e1[X] - n[X]*e1[Z];
    e2[Z] = n[X]*e1[Y] - n[Y]*e1[X];

    // Inclusion center, PBC-adjusted near planePt (or pin)
    real8 ctr[3] = { inc->position[X], inc->position[Y], inc->position[Z] };
    real8 ref[3] = { planePt[X], planePt[Y], planePt[Z] };
    PBCPOSITION(param, ref[X], ref[Y], ref[Z], &ctr[X], &ctr[Y], &ctr[Z]);

    real8 aR = inc->radius[X];
    real8 bR = inc->radius[Y];
    real8 cR = inc->radius[Z];
    if ((aR <= 0.0) || (bR <= 0.0) || (cR <= 0.0)) return 0;

    // Build ellipsoid quadratic form A in global coords:
    // A = u1u1^T/a^2 + u2u2^T/b^2 + u3u3^T/c^2
    real8 u1[3], u2[3], u3[3];
    u1[X] = inc->rotation[0][0]; u1[Y] = inc->rotation[0][1]; u1[Z] = inc->rotation[0][2];
    u2[X] = inc->rotation[1][0]; u2[Y] = inc->rotation[1][1]; u2[Z] = inc->rotation[1][2];
    u3[X] = u1[Y]*u2[Z] - u1[Z]*u2[Y];
    u3[Y] = u1[Z]*u2[X] - u1[X]*u2[Z];
    u3[Z] = u1[X]*u2[Y] - u1[Y]*u2[X];

    real8 ia2 = 1.0/(aR*aR);
    real8 ib2 = 1.0/(bR*bR);
    real8 ic2 = 1.0/(cR*cR);

    real8 A[3][3];
    for (int i=0;i<3;i++) for (int j=0;j<3;j++) {
        real8 u1i = (i==0)?u1[X]:(i==1)?u1[Y]:u1[Z];
        real8 u1j = (j==0)?u1[X]:(j==1)?u1[Y]:u1[Z];
        real8 u2i = (i==0)?u2[X]:(i==1)?u2[Y]:u2[Z];
        real8 u2j = (j==0)?u2[X]:(j==1)?u2[Y]:u2[Z];
        real8 u3i = (i==0)?u3[X]:(i==1)?u3[Y]:u3[Z];
        real8 u3j = (j==0)?u3[X]:(j==1)?u3[Y]:u3[Z];
        A[i][j] = ia2*u1i*u1j + ib2*u2i*u2j + ic2*u3i*u3j;
    }

    // Ensure pin lies on plane (project once, just in case)
    real8 pinP[3] = { pin[X], pin[Y], pin[Z] };
    real8 dp = (pinP[X]-planePt[X])*n[X] + (pinP[Y]-planePt[Y])*n[Y] + (pinP[Z]-planePt[Z])*n[Z];
    pinP[X] -= dp*n[X]; pinP[Y] -= dp*n[Y]; pinP[Z] -= dp*n[Z];

    // 2D coordinates of pin in (e1,e2) basis relative to planePt
    real8 r0[3] = { pinP[X]-planePt[X], pinP[Y]-planePt[Y], pinP[Z]-planePt[Z] };
    real8 u0 = r0[X]*e1[X] + r0[Y]*e1[Y] + r0[Z]*e1[Z];
    real8 v0 = r0[X]*e2[X] + r0[Y]*e2[Y] + r0[Z]*e2[Z];

    // Build quadratic F(u,v) = d^T A d - 1 = 0 where d = (planePt-ctr) + u e1 + v e2
    real8 d0[3] = { planePt[X]-ctr[X], planePt[Y]-ctr[Y], planePt[Z]-ctr[Z] };

    // Helper dot products with A
    auto Ax = [&](const real8 x[3], real8 out[3]){
        out[X] = A[0][0]*x[X] + A[0][1]*x[Y] + A[0][2]*x[Z];
        out[Y] = A[1][0]*x[X] + A[1][1]*x[Y] + A[1][2]*x[Z];
        out[Z] = A[2][0]*x[X] + A[2][1]*x[Y] + A[2][2]*x[Z];
    };

    real8 Ae1[3], Ae2[3], Ad0[3];
    Ax(e1, Ae1); Ax(e2, Ae2); Ax(d0, Ad0);

    // Q (2x2), q (2), s (scalar)
    real8 Q11 = e1[X]*Ae1[X] + e1[Y]*Ae1[Y] + e1[Z]*Ae1[Z];
    real8 Q12 = e1[X]*Ae2[X] + e1[Y]*Ae2[Y] + e1[Z]*Ae2[Z];
    real8 Q22 = e2[X]*Ae2[X] + e2[Y]*Ae2[Y] + e2[Z]*Ae2[Z];
    real8 q1  = e1[X]*Ad0[X] + e1[Y]*Ad0[Y] + e1[Z]*Ad0[Z];
    real8 q2  = e2[X]*Ad0[X] + e2[Y]*Ad0[Y] + e2[Z]*Ad0[Z];
    real8 s0  = d0[X]*Ad0[X] + d0[Y]*Ad0[Y] + d0[Z]*Ad0[Z] - 1.0;

    auto F = [&](real8 u, real8 v)->real8 {
        return Q11*u*u + 2.0*Q12*u*v + Q22*v*v + 2.0*q1*u + 2.0*q2*v + s0;
    };

    real8 F0 = F(u0,v0);
    if (fabs(F0) < 1.0e-12) {
        pout[X] = pinP[X]; pout[Y] = pinP[Y]; pout[Z] = pinP[Z];
        return 1;
    }

    // Solve KKT via scalar lambda with robust bisection:
    // (I + λQ) z = z0 - λq, and enforce F(z(λ)) = 0
    auto Zof = [&](real8 lam, real8 &u, real8 &v)->int {
        real8 M11 = 1.0 + lam*Q11;
        real8 M12 =       lam*Q12;
        real8 M22 = 1.0 + lam*Q22;
        real8 det = M11*M22 - M12*M12;
        if (fabs(det) < 1.0e-20) return 0;

        real8 b1 = u0 - lam*q1;
        real8 b2 = v0 - lam*q2;

        u = ( M22*b1 - M12*b2)/det;
        v = (-M12*b1 + M11*b2)/det;
        return 1;
    };

    auto G = [&](real8 lam)->real8 {
        real8 u,v;
        if (!Zof(lam,u,v)) return (F0>0? 1.0e30 : -1.0e30);
        return F(u,v);
    };

    // Bracket root
    real8 lo, hi;
    if (F0 > 0.0) { // outside -> need positive lambda
        lo = 0.0; hi = 1.0;
        for (int k=0;k<80;k++){
            real8 Gh = G(hi);
            if (Gh < 0.0) break;
            hi *= 2.0;
            if (hi > 1.0e12) return 0;
        }
    } else { // inside -> need negative lambda
        hi = 0.0; lo = -1.0;
        for (int k=0;k<80;k++){
            real8 Gl = G(lo);
            if (Gl > 0.0) break;
            lo *= 2.0;
            if (lo < -1.0e12) return 0;
        }
    }

    // Bisection
    real8 lam = 0.0;
    for (int it=0; it<80; it++){
        lam = 0.5*(lo+hi);
        real8 Gm = G(lam);
        if (fabs(Gm) < 1.0e-12) break;
        if (Gm > 0.0) lo = lam; else hi = lam;
    }

    // Final point
    real8 uf,vf;
    if (!Zof(lam, uf, vf)) return 0;

    // Put on plane
    pout[X] = planePt[X] + uf*e1[X] + vf*e2[X];
    pout[Y] = planePt[Y] + uf*e1[Y] + vf*e2[Y];
    pout[Z] = planePt[Z] + uf*e1[Z] + vf*e2[Z];

    // Slightly push outward along in-plane gradient to avoid re-penetration
    // (optional, very small)
    real8 d[3] = { pout[X]-ctr[X], pout[Y]-ctr[Y], pout[Z]-ctr[Z] };
    real8 Ad[3]; Ax(d, Ad);
    // in-plane component of grad
    real8 gm = Ad[X]*n[X] + Ad[Y]*n[Y] + Ad[Z]*n[Z];
    Ad[X] -= gm*n[X]; Ad[Y] -= gm*n[Y]; Ad[Z] -= gm*n[Z];
    real8 gMag = sqrt(Ad[X]*Ad[X] + Ad[Y]*Ad[Y] + Ad[Z]*Ad[Z]);
    if (gMag > 1.0e-20) {
        real8 eps = 1.0e-6;
        pout[X] += eps*(Ad[X]/gMag);
        pout[Y] += eps*(Ad[Y]/gMag);
        pout[Z] += eps*(Ad[Z]/gMag);
    }

    FoldBox(param, &pout[X], &pout[Y], &pout[Z]);
    return 1;
}
/*---------------------------------------------------------------------------
 *
 *      Function:     IntersectionSphereSegment
 *      Description:  Determine the number of intersecting points between
 *                    a sphere and a line segment, the fraction of the
 *                    segment inside the sphere, the the derivative of
 *                    these fractions.
 *
 *      Note : Before entering this routine, care must be made to ensure 
 *             that PBC is taken into account for pos1, pos2, and C before 
 *             rotation in the case of an ellipse.
 * 
 *      Arguments:
 *          p1 : node starting the segment
 *          p2 : node ending the segment
 *          C  : sphere center
 *          r  : sphere radius
 *          n  : number of intersection points
 *          ratio : fractions of the line segments where the intersections points
 *                   are. These fractions can be outside the range of [0, 1].
 *
 *-------------------------------------------------------------------------*/
static int IntersectionSphereSegment(real8 p1[3], real8 p2[3],real8 C[3], 
                                     real8 r, int *n, real8 ratio[2])
{
   real8 v[3], t[3], R[2][3], s[2], nd[3], sinter[2];  
   real8 invL, d2, test, sdiff;
   (*n) = 0;

   int isInter = 0;
   ratio[0] = 0.0;
   ratio[1] = 0.0;

   v[X] = p2[X]-p1[X];
   v[Y] = p2[Y]-p1[Y];
   v[Z] = p2[Z]-p1[Z];

   invL = 1./Normal(v);
   t[X] = v[X]*invL;
   t[Y] = v[Y]*invL;
   t[Z] = v[Z]*invL;

   R[0][X] = p1[X] - C[X];
   R[0][Y] = p1[Y] - C[Y];
   R[0][Z] = p1[Z] - C[Z];

   R[1][X] = p2[X] - C[X];
   R[1][Y] = p2[Y] - C[Y];
   R[1][Z] = p2[Z] - C[Z];

   s[0] = DotProduct(R[0],t);
   s[1] = DotProduct(R[1],t);

   nd[X] = R[0][X] - s[0]*t[X];
   nd[Y] = R[0][Y] - s[0]*t[Y];
   nd[Z] = R[0][Z] - s[0]*t[Z];

   d2 = DotProduct(nd,nd);
   test = r*r-d2;

   if (test < 0.0) return 0;
   
   sinter[0] = -sqrt(test);
   sinter[1] = -sinter[0];
   
   sdiff = s[1]-s[0];
   ratio[0] = (sinter[0] - s[0])/sdiff;
   ratio[1] = (sinter[1] - s[0])/sdiff;

   if (ratio[0]>0.0 && ratio[0]<1) (*n)++;
   if (ratio[1]>0.0 && ratio[1]<1) (*n)++;
   if ( (*n) > 0) 
      isInter = 1;
   else
      isInter = 0;

   // The whole segment could be inside the particle.
   // This counts as an intersection.
   if (ratio[0]<0.0 && ratio[1]>1) isInter = 1;

   ratio[0] = MIN(1.0, MAX(ratio[0], 0.0));
   ratio[1] = MIN(1.0, MAX(ratio[1], 0.0));
   
   return isInter;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     TransformEllips2Sphere
 *      Description:  Matrix defining the rotation and scaling of coordinate
 *                    system from an ellipsoid to a sphere. Used in the 
 *                    mobility laws called Eshelby for now.
 *                    The matrix B is defined as the product between the diagonal 
 *                    matrix of the semi-principal axis multiplied with the rotation
 *                    matrix of these axis
 * 
 *      Arguments:
 *            rotation[2][3] : two vectors defining the rotation matrix
 *            radius[3]      : semi-principal axis
 *             
 *            B[3][3]        : The matrix defining the transformation from ellipoidal
 *                             coordinates to spherical coordinates.
 *
 *---------------------------------------------------------------------------*/
static void TransformEllips2Sphere(real8 rotation[2][3], real8 radius[3], real8 B[3][3])
{
  real8 vec[3];
  NormalizedCrossVector(rotation[0],rotation[1],vec);

  B[0][0] = rotation[0][0]/radius[0];
  B[0][1] = rotation[0][1]/radius[0];
  B[0][2] = rotation[0][2]/radius[0];
  
  B[1][0] = rotation[1][0]/radius[1];
  B[1][1] = rotation[1][1]/radius[1];
  B[1][2] = rotation[1][2]/radius[1];
  
  B[2][0] = vec[0]/radius[2];
  B[2][1] = vec[1]/radius[2];
  B[2][2] = vec[2]/radius[2];

}

/*---------------------------------------------------------------------------
 *
 * Function:     GetMinDistPointToEllipsoid
 * Description:  Calculates the minimum squared distance from a point P
 * to the surface of an ellipsoid using Newton's Method.
 * Based on Literature Chapter 5, Eq (5.20) - (5.22).
 *
 * Arguments:
 * p:        Query point (global coordinates)
 * c:        Ellipsoid center
 * radius:   Ellipsoid semi-axes (r1, r2, r3)
 * Rotation: Rotation matrix vectors (u1, u2)
 *
 * Returns:      Minimum distance squared
 *
 *-------------------------------------------------------------------------*/
static real8 GetMinDistPointToEllipsoid(real8 p[3], real8 c[3], 
                                        real8 radius[3], real8 Rotation[2][3])
{
    // 1. 构建完整的局部坐标系旋转矩阵 (3x3)
    // Rotation 只提供了两个轴 u1, u2，我们需要计算 u3 = u1 x u2
    real8 u1[3], u2[3], u3[3];
    VECTOR_COPY(u1, Rotation[0]);
    VECTOR_COPY(u2, Rotation[1]);
    NormalizedCrossVector(u1, u2, u3); 

    // 2. 将点 P 转换到椭球的局部主轴坐标系 (Local Frame)
    // P_local = Rot^T * (P - C)
    real8 diff[3], p_local[3];
    diff[X] = p[X] - c[X];
    diff[Y] = p[Y] - c[Y];
    diff[Z] = p[Z] - c[Z];

    // 投影到三个主轴上 (相当于乘以旋转矩阵的转置)
    p_local[X] = DotProduct(diff, u1);
    p_local[Y] = DotProduct(diff, u2);
    p_local[Z] = DotProduct(diff, u3);

    // 3. 牛顿迭代求解最近点 (Newton's Method)
    // 求解方程组 F(y, lambda) = 0
    // 变量 x = [y1, y2, y3, lambda]
    
    // 初始猜测: 将局部点投影到球面上 (文献 Eq 5.22 的简化版)
    real8 y[3];
    real8 r2[3] = {radius[0]*radius[0], radius[1]*radius[1], radius[2]*radius[2]};
    real8 r_inv2[3] = {1.0/r2[0], 1.0/r2[1], 1.0/r2[2]};
    
    // 简单初始化：直接取 p_local
    VECTOR_COPY(y, p_local);
    
    real8 lambda = 0.0;
    const int MAX_ITER = 10;
    const real8 TOL = 1.0e-6; // 收敛容差

    for (int iter = 0; iter < MAX_ITER; iter++) {
        // 构建残差向量 F (4x1)
        // F_i = y_i/r_i^2 + 2*y_i*lambda - p_local_i/r_i^2 (文献公式可能有符号差异，这里推导修正)
        // 拉格朗日函数 L = |y-p|^2 + lambda * (y^2/r^2 - 1)
        // dL/dy = 2(y-p) + 2*lambda*y/r^2 = 0  => y - p + lambda*y/r^2 = 0
        // y * (1 + lambda/r^2) = p  => y = p / (1 + lambda/r^2)
        // 这种形式更适合直接迭代 lambda，但为了严格遵循牛顿法求解非线性方程组：
        
        real8 F[4];
        for(int i=0; i<3; i++) {
            // 方程: y_i + y_i * lambda / r_i^2 - p_local_i = 0  (由梯度导数简化)
            // 这里 lambda 吸收了常数因子2
            F[i] = y[i] * (1.0 + lambda * r_inv2[i]) - p_local[i];
        }
        F[3] = (y[0]*y[0]*r_inv2[0] + y[1]*y[1]*r_inv2[1] + y[2]*y[2]*r_inv2[2]) - 1.0;

        // 检查收敛
        if (fabs(F[0])<TOL && fabs(F[1])<TOL && fabs(F[2])<TOL && fabs(F[3])<TOL) break;

        // 构建 Jacobian 矩阵 J (4x4)
        real8 J[4][4];
        Matrix44_Zero(J);
        
        for(int i=0; i<3; i++) {
            // d(F_i) / d(y_i) = 1 + lambda / r_i^2
            J[i][i] = 1.0 + lambda * r_inv2[i];
            // d(F_i) / d(lambda) = y_i / r_i^2
            J[i][3] = y[i] * r_inv2[i];
            // d(F_3) / d(y_i) = 2 * y_i / r_i^2
            J[3][i] = 2.0 * y[i] * r_inv2[i];
        }
        
        // 求解 J * dX = -F
        real8 invJ[4][4];
        if (Matrix44_Inverse(invJ, J) == -1) {
            // 如果矩阵奇异（极少见），稍微扰动一下继续
            lambda += 0.1;
            continue;
        }
        
        real8 dX[4] = {0,0,0,0};
        // dX = invJ * (-F)
        for(int i=0; i<4; i++) {
            for(int j=0; j<4; j++) {
                dX[i] -= invJ[i][j] * F[j];
            }
        }
        
        // 更新状态
        y[0] += dX[0];
        y[1] += dX[1];
        y[2] += dX[2];
        lambda += dX[3];
    }

    // 计算距离平方: |p_local - y_surface|^2
    real8 dist2 = (p_local[0]-y[0])*(p_local[0]-y[0]) + 
                  (p_local[1]-y[1])*(p_local[1]-y[1]) + 
                  (p_local[2]-y[2])*(p_local[2]-y[2]);
                  
    return dist2;
}
/*---------------------------------------------------------------------------
 *
 *      Function:     IntersectionEllipseSegment
 *      Description:  Determine the number of intersecting points between the sphere
 *                    and the line segment and the fraction of the segment inside the 
 *                    sphere.
 *
 *      Note : Before entering this routine, care must be made to ensure that PBC is 
 *             taken into account for pos1, pos2, and C before rotation in the case 
 *             of an ellipse.
 * 
 *      Arguments:
 *          C         : position of the ellipse
 *          Radius    : Sphere radius if ellispe is a sphere. Set to 1 in the case of 
 *                      an ellipse
 *          Rotation  : Two vectors defining the rotation matrix
 *          p1        : node starting the segment
 *          p2        : node ending the segment
 *          ninc      : number of intersection points
 *          ratio     : fractions of the line segments where the intersections points
 *                      are. These fractions can be outside the range of [0, 1].
 *
 *-------------------------------------------------------------------------*/
int IntersectionEllipseSegment(real8  C[3], real8 Radius[3], real8 Rotation[2][3],
                               real8 p1[3], real8 p2[3], real8 r, int *ninc, 
                               real8 ratio[2])
{
   int IsIntersect = 0;
   real8 B[3][3], Rp1[3], Rp2[3], RC[3];

   // Transform the ellipse into a sphere to get fractions of dislocation line segment
   // inside the ellipse
   TransformEllips2Sphere(Rotation, Radius, B);
   Matrix33Vector3Multiply(B,p1,Rp1);
   Matrix33Vector3Multiply(B,p2,Rp2);
   Matrix33Vector3Multiply(B,C,RC);
   IsIntersect = IntersectionSphereSegment(Rp1, Rp2, RC, r, ninc, ratio);
   return IsIntersect;
}

/*---------------------------------------------------------------------------
 *
 * Function:     IntersectionInformation
 * Description:  Modified to include precise Ellipsoid proximity check.
 *
 *-------------------------------------------------------------------------*/
int IntersectionInformation(Home_t *home, real8 pos1[3], real8 pos2[3], EInclusion_t *inclusion, 
                            real8 newPos1[3], real8 newPos2[3], int *ninc, real8 ratio[2])
{
   int IsIntercept=0;
   real8 cellCtr[3];

   Param_t *param;
   Cell_t *cell;

   param = home->param;

   int cellID;
   real8 C[3], radius[3], Rotation[2][3];
   cellID = inclusion->cellID;
   VECTOR_COPY(C,inclusion->position);
   VECTOR_COPY(radius,inclusion->radius);
   VECTOR_COPY(Rotation[0],inclusion->rotation[0]);
   VECTOR_COPY(Rotation[1],inclusion->rotation[1]);
  
   // Ensure PBC (Period Boundary Conditions)
   newPos1[X] = pos1[X]; newPos1[Y] = pos1[Y]; newPos1[Z] = pos1[Z];
   newPos2[X] = pos2[X]; newPos2[Y] = pos2[Y]; newPos2[Z] = pos2[Z];
   
   cell = LookupCell(home, cellID);
   
   if (cell == (Cell_t *)NULL) {
      FindCellCenter(param, C[X], C[Y], C[Z], 1,
                     &cellCtr[X], &cellCtr[Y], &cellCtr[Z]);
   } else {
      cellCtr[X] = cell->center[X];
      cellCtr[Y] = cell->center[Y];
      cellCtr[Z] = cell->center[Z];
   }
   
   PBCPOSITION(param, cellCtr[X], cellCtr[Y], cellCtr[Z],
               &newPos1[X], &newPos1[Y], &newPos1[Z]);
               
    // 2. 【新增】修正 C (析出相中心)，让它也靠近 cellCtr
   PBCPOSITION(param, cellCtr[X], cellCtr[Y], cellCtr[Z],
               &C[X], &C[Y], &C[Z]);

   PBCPOSITION(param, newPos1[X], newPos1[Y], newPos1[Z],
               &newPos2[X], &newPos2[Y], &newPos2[Z]);
   
   // 1. 首先尝试原有的“穿透”检测 (Transformation Method)
   // 这个方法对于判断 "Inside" 很有效
   IsIntercept = IntersectionEllipseSegment(C, radius, Rotation, newPos1, newPos2, 
                                            1.0, ninc, ratio);

   // 2. 如果未检测到穿透 (IsIntercept == 0)，则使用牛顿法检测“贴面接触” (Surface Contact)
   // 这是为了实现 Orowan 绕过机制：位错线在外部紧贴表面
   if (IsIntercept == 0) {
       // 定义接触容差 (Contact Tolerance)
       // 通常取 0.5b 到 1.0b，这取决于你的网格精度
       real8 contactTol = 0.5; // units: b (simulation length unit); contact tolerance ~0.5b
       real8 contactTol2 = contactTol * contactTol;

       // 检查端点 1 的距离
       real8 dist2_p1 = GetMinDistPointToEllipsoid(newPos1, C, radius, Rotation);
       
       // 检查端点 2 的距离
       real8 dist2_p2 = GetMinDistPointToEllipsoid(newPos2, C, radius, Rotation);

       // 如果任意一个端点距离表面足够近，我们强制认为发生了碰撞
       if (dist2_p1 < contactTol2 || dist2_p2 < contactTol2) {
           IsIntercept = 1;
           // 这里我们不修改 ratio，因为不是穿透，而是接触
           // 这会触发 HandleParticleCollisions 里的逻辑
       }
   }

   return IsIntercept;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     GetIndexFromInclusionID()
 *      Description:  Locate the specified inclusion in the local domain's
 *                    inclusion list and return the index in the list
 *                    to the caller.  
 *
 *                    NOTE: The inclusions are NOT sorted, so unless we
 *                    modify the code to do so, we have to do a linear
 *                    search through the list for the inclusion.
 *
 *      Arguments:
 *          inclusionID  self-explanatory
 *
 *      Returns:  The index (0-based) in the local inclusion list of the
 *                inclusion with the provided ID.
 *
 *-------------------------------------------------------------------------*/
int GetIndexFromInclusionID(Home_t *home, int inclusionID)
{
        int                numEntries;
        SegPartIntersect_t *entry, *tmpList;

        // Access to segment/particle intersection list
        numEntries = home->segpartIntersectCnt;
        tmpList = home->segpartIntersectList;

        // Loop over all entries in the seg/part list
        for (int i = 0; i <numEntries ; i++) 
        {
           entry = &tmpList[i];           
           if (entry == (SegPartIntersect_t *)NULL)  return(-1);
           
           for (int j = 0; j < (entry)->numIntersections ; j++) 
           {
              if (inclusionID == entry->inclusionID[j])
              {
                 return (entry->inclusionIndex[j]);
              }
           }
        }

        return(-1);
}


/*---------------------------------------------------------------------------
 *
 *      Function:     SegPartListAppend()
 *      Description:  Append the provided segment/particle intersection
 *                    information to the segment/particle intersection
 *                    list.  Caller must resort the list before using
 *                    the SegPartListLookup() function.
 *
 *      Arguments:
 *          newInfo  Pointer to structure particle intersection data
 *                   for a segment.  Data is copied from this structure
 *                   into the segment/particle intersction list.
 *
 *-------------------------------------------------------------------------*/
void SegPartListAppend(Home_t *home, SegPartIntersect_t *newInfo)
{
        int                i, numEntries, allocedEntries;
        SegPartIntersect_t *intersectInfo, *tmpList;

/*
 *      The segment/particle intersection list is globally
 *      available to all threads, so we have to ensure only
 *      one thread accesses it at a time
 */
#ifdef _OPENMP
#pragma omp critical (ACCESS_SEG_PART_LIST)
#endif
        {
            numEntries = home->segpartIntersectCnt;
            allocedEntries = home->segpartIntersectListSize;
            tmpList = home->segpartIntersectList;

/*
 *          Make sure there is enough space on the list for a segment
 */
            if (numEntries >= allocedEntries) {
                allocedEntries += 100;
                tmpList = (SegPartIntersect_t *)realloc(tmpList,
                           allocedEntries * sizeof(SegPartIntersect_t));
            }

/*
 *          Get a pointer to the next free entry in the list and
 *          initialize it.  IMPORTANT: Make sure that the node tags
 *          identifying the segment are stored such that the first tag
 *          for a segment is "less than" the second as determined
 *          by OrderTags().  This is so we're consistent in the manner
 *          we sort the list and search for entries on the list later.
 */
            intersectInfo = &tmpList[numEntries];
            memset(intersectInfo, 0, sizeof(SegPartIntersect_t));

            if (OrderTags(&newInfo->tag1, &newInfo->tag2) < 0) {
                intersectInfo->tag1.domainID = newInfo->tag1.domainID;
                intersectInfo->tag1.index    = newInfo->tag1.index;
                intersectInfo->tag2.domainID = newInfo->tag2.domainID;
                intersectInfo->tag2.index    = newInfo->tag2.index;
            } else {
                intersectInfo->tag1.domainID = newInfo->tag2.domainID;
                intersectInfo->tag1.index    = newInfo->tag2.index;
                intersectInfo->tag2.domainID = newInfo->tag1.domainID;
                intersectInfo->tag2.index    = newInfo->tag1.index;
            }

/*
 *          Probably don't need to copy the inclusion IDs here since
 *          this function should only be called when reading the intersect
 *          data from a remote domain and will not be sending it any
 *          further, but...
 */
            for (i = 0; i < newInfo->numIntersections; i++) {
                intersectInfo->inclusionIndex[i] = newInfo->inclusionIndex[i];
                intersectInfo->inclusionID[i] = newInfo->inclusionID[i];
            }

            home->segpartIntersectCnt = numEntries + 1;
            home->segpartIntersectListSize = allocedEntries;
            home->segpartIntersectList = tmpList;

        }  /* end "omp critical" section */

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     FindSegPartIntersects
 *      Description:  Determine which inclusions (particles) intersect
 *                    the specified segment and add the info to the
 *                    seg/particle intersection list.  
 * 
 *                    NOTE:  This function does NOT verify that the
 *                           list does not already contain the intersection
 *                           data for the specified segment.  It is up to
 *                           the caller to ensure he does not call this
 *                           function for the same segment multiple times.
 *
 *      Arguments:
 *          tag1    Tag of one of the segment endpoints
 *          tag2    Tag for the other segment endpoint
 *          pos1    Coordinates of the node associated with tag1
 *          pos2    Coordinates of the node associated with tag2, adjusted
 *                  to the periodic image closest to pos1.
 *          cellID  ID of the cell containing pos1
 *          cellCtr Coordinates of the center of cell <cellID>
 *
 *-------------------------------------------------------------------------*/
static void FindSegPartIntersects(Home_t *home, Tag_t *tag1, Tag_t *tag2,
                                  real8 pos1[3], real8 pos2[3], int cellID,
                                  real8 cellCtr[3])
{
        int                j, numNbrCells, inclusionIndex, inclusionCellID;
        real8              sinter[2];
        real8              incPos[3];
        Cell_t             *cell, *inclusionCell;
        Param_t            *param;
        EInclusion_t       *inclusion;
        SegPartIntersect_t *intersection = (SegPartIntersect_t *)NULL;

        param = home->param;
        cell = LookupCell(home, cellID);
        numNbrCells = cell->nbrCount;
        intersection = (SegPartIntersect_t *)NULL;

/*
 *      First loop over all neighboring cells, then on last loop
 *      iteration do the current cell.
 */
        for (j = 0; j <= numNbrCells; j++) {

            if (j < numNbrCells) {
                inclusionCellID = cell->nbrList[j];
            } else {
                inclusionCellID = cellID;
            }

/*
 *          It *may* be possible that we don't have cell info for
 *          all the cells neighboring the provided <cellID>.  In
 *          that case we just do the best we can.
 */
            inclusionCell = LookupCell(home, inclusionCellID);
            if (inclusionCell == (Cell_t *)NULL) continue;

/*
 *          NOTE: If the cell is a ghost cell, alter the cellID and cell
 *          pointer to the corresponding base cell.  Needed in order to pick
 *          up the proper inclusionQ.
 */
            inclusionCell = LookupCell(home, inclusionCellID);

            if (inclusionCell->baseIdx >= 0) {
                inclusionCellID = inclusionCell->baseIdx;
                inclusionCell = LookupCell(home, inclusionCellID);
                if (inclusionCell == (Cell_t *)NULL) continue;
            }

            inclusionIndex = inclusionCell->inclusionQ;

/*
 *          Loop through all the inclusion in this cell
 */
            while (inclusionIndex >= 0) {

                int ninc;

                inclusion = &home->eshelbyInclusions[inclusionIndex];


/*
 *              Adjust the center of the inclusion to that of the periodic
 *              image closest to the center of the cell owning the first node.
 */
                incPos[X] = inclusion->position[X];
                incPos[Y] = inclusion->position[Y];
                incPos[Z] = inclusion->position[Z];

                PBCPOSITION(param, cellCtr[X], cellCtr[Y], cellCtr[Z],
                            &incPos[X], &incPos[Y], &incPos[Z]);

                int IsIntercept = 0;
                IsIntercept = IntersectionEllipseSegment(incPos, inclusion->radius, inclusion->rotation,
                                                         pos1, pos2, 1.0, &ninc, sinter);

/*
 *              If the segment does not intersect the inclusion, i.e number 
 *              of intersecting points is zero, skip it
 */
                if (IsIntercept)
                {
                   SegPartListUpdate(home, &intersection, tag1, tag2,
                                     inclusionIndex, inclusion->id);
                }

                inclusionIndex = inclusion->nextInCell;
            }

        }  /* loop over cells */

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     FindNodePartIntersects
 *      Description:  Given a node, find out which inclusions are intersected
 *                    by any of the segments attached to the node.  The
 *                    segment/particle intersect list will be updated with
 *                    the info.
 *
 *                    NOTE:  It is up to the caller to sort the seg/particle
 *                    intersection list after he's done generating it.
 *
 *
 *      Returns:   Pointer to structure of segment/particle intersection
 *                 information, or NULL if segment does not intersect any
 *                 inclusions.
 *
 *-------------------------------------------------------------------------*/
void FindNodePartIntersects(Home_t *home, Node_t *node)
{
        int     i, cellID;
        real8   eps = 1.0e-6, dx, dy, dz;
        real8   pos1[3], pos2[3], cellCtr[3];
        Node_t  *nbrNode;
        Cell_t  *cell;
        SegPartIntersect_t *intersection;

        pos1[X] = node->x;
        pos1[Y] = node->y;
        pos1[Z] = node->z;

        cellID = node->cellIdx;
        cell = LookupCell(home, cellID);

        cellCtr[X] = cell->center[X];
        cellCtr[Y] = cell->center[Y];
        cellCtr[Z] = cell->center[Z];

        for (i = 0; i < node->numNbrs; i++) {

            if ((nbrNode = GetNeighborNode(home, node, i)) == (Node_t *)NULL) {
                continue;
            }

            dx = nbrNode->x - node->x;
            dy = nbrNode->y - node->y;
            dz = nbrNode->z - node->z;
            
            ZImage(home->param, &dx, &dy, &dz);


            pos2[X] = node->x + dx;
            pos2[Y] = node->y + dy;
            pos2[Z] = node->z + dz;


/*
 *          Zero length arms can occur when we're testing multinode splits.
 *          Ignore such segments.
 */
            if ((fabs(pos1[X]-pos2[X]) < eps) &&
                (fabs(pos1[Y]-pos2[Y]) < eps) &&
                (fabs(pos1[Z]-pos2[Z]) < eps)) {
                continue;
            }

/*
 *          If this function is called for two connected nodes, we have to
 *          be sure we don't search for intersecting particles multiple
 *          times, so only add the segment/particle intersections if the
 *          segment is not yet on the list.  (Note: list may be unsorted at
 *          this point, so we have to use a lookup function that handles
 *          an unsorted list)
 */
            intersection = SegPartListUnsortedLookup(home, &node->myTag,
                                                     &nbrNode->myTag);
            if (intersection != (SegPartIntersect_t *)NULL) {
                continue;
            }

            FindSegPartIntersects(home, &node->myTag, &nbrNode->myTag,
                                  pos1, pos2, cellID, cellCtr);
            
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     SegPartListCmp()
 *      Description:  Compares two SegPartIntersect_t structures based
 *                    on the <tag1> and <tag2> values.  This function
 *                    is compatible with use by system sorting and
 *                    searching functions such as qsort(), bsearch(), etc.
 *
 *      Returns:        -1 if  a <  b
 *                       0 if  a == b
 *                       1 if  a >  b
 *
 *-------------------------------------------------------------------------*/
static int SegPartListCmp(const void *a, const void *b)
{
        SegPartIntersect_t *info1 = (SegPartIntersect_t *)a;
        SegPartIntersect_t *info2 = (SegPartIntersect_t *)b;

/*
 *      First compare first tag 
 */
        if (info1->tag1.domainID < info2->tag1.domainID) return(-1);
        if (info1->tag1.domainID > info2->tag1.domainID) return(1);

        if (info1->tag1.index < info2->tag1.index) return(-1);
        if (info1->tag1.index > info2->tag1.index) return(1);

/*
 *      Both entries have the same first tag, so compare the second
 *      tag of each.
 */
        if (info1->tag2.domainID < info2->tag2.domainID) return(-1);
        if (info1->tag2.domainID > info2->tag2.domainID) return(1);

        if (info1->tag2.index < info2->tag2.index) return(-1);
        if (info1->tag2.index > info2->tag2.index) return(1);

        return(0);
}


/*---------------------------------------------------------------------------
 *
 *      Function:     SegPartListSort()
 *      Description:  Sort the list of particle-intersecting segments
 *
 *-------------------------------------------------------------------------*/
void SegPartListSort(Home_t *home)
{
        if (home->segpartIntersectCnt > 1) {
            qsort(home->segpartIntersectList, home->segpartIntersectCnt,
                  sizeof(SegPartIntersect_t), SegPartListCmp);
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     SegPartListLookup()
 *      Description:  Given a pair of node tags defining a segment, return
 *                    to the caller the corresponding entry (if any) in the
 *                    segment/particle intersection list.
 *
 *      Arguments:
 *          tag1  Node tag of one of the segment endpoints
 *          tag2  Node tag of the other segment endpoint
 *
 *      Returns:  Pointer to the entry in the segment/particle
 *                intersection list with endpoints <tag1> and <tag2>
 *                or NULL if the segment is not found on the list
 *
 *-------------------------------------------------------------------------*/
SegPartIntersect_t *SegPartListLookup(Home_t *home, Tag_t *tag1, Tag_t *tag2)
{
        SegPartIntersect_t key;
        SegPartIntersect_t *entry = (SegPartIntersect_t *)NULL;

/*
 *      Node tags in the intersection list are stored such that the
 *      first tag for a segment is "less than" the second as determined
 *      by OrderTags(), so make sure we set up the tags for the search
 *      key the same way.
 */
        if (OrderTags(tag1, tag2) < 0) {
            key.tag1.domainID = tag1->domainID;
            key.tag1.index    = tag1->index;
            key.tag2.domainID = tag2->domainID;
            key.tag2.index    = tag2->index;
        } else {
            key.tag1.domainID = tag2->domainID;
            key.tag1.index    = tag2->index;
            key.tag2.domainID = tag1->domainID;
            key.tag2.index    = tag1->index;
        }

        entry = (SegPartIntersect_t *)bsearch(&key, home->segpartIntersectList,
                                              home->segpartIntersectCnt,
                                              sizeof(SegPartIntersect_t),
                                              SegPartListCmp);

        return(entry);
}


/*---------------------------------------------------------------------------
 *
 *      Function:     SegPartListUnsortedLookup()
 *      Description:  Given a pair of node tags defining a segment, return
 *                    to the caller the corresponding entry (if any) in the
 *                    segment/particle intersection list.  This function
 *                    does NOT assume the list is currently sorted and 
 *                    does a linear search through the array.  This is
 *                    only intended to be used when dealing with a very
 *                    limited set of segments, and so shouldn't be very
 *                    expensive.
 *
 *      Arguments:
 *          tag1  Node tag of one of the segment endpoints
 *          tag2  Node tag of the other segment endpoint
 *
 *      Returns:  Pointer to the entry in the segment/particle
 *                intersection list with endpoints <tag1> and <tag2>
 *                or NULL if the segment is not found on the list
 *
 *-------------------------------------------------------------------------*/
SegPartIntersect_t *SegPartListUnsortedLookup(Home_t *home, Tag_t *tag1,
                                              Tag_t *tag2)
{
        int                i;
        Tag_t              *tmpTag1, *tmpTag2;
        SegPartIntersect_t *entry = (SegPartIntersect_t *)NULL;

        if (OrderTags(tag1, tag2) < 0) {
           tmpTag1 = tag1;
           tmpTag2 = tag2;
        } else {
           tmpTag2 = tag1;
           tmpTag1 = tag2;
        }

        for (i = 0; i < home->segpartIntersectCnt; i++) {

            entry = &home->segpartIntersectList[i];

            if ((tmpTag1->domainID == entry->tag1.domainID) &&
                (tmpTag1->index    == entry->tag1.index   ) &&
                (tmpTag2->domainID == entry->tag2.domainID) &&
                (tmpTag2->index    == entry->tag1.index   )) {

                return(entry);
            }
        }

        return((SegPartIntersect_t *)NULL);
}


/*---------------------------------------------------------------------------
 *
 *      Function:     SegPartListUpdate()
 *      Description:  Allocate (if necessary) a new entry in the segment/
 *                    particle intersection list, initialize the
 *                    segment information, and return a pointer to the
 *                    entry to the caller.
 *
 *      Arguments:
 *          entry          Location containing pointer to the intersection
 *                         data to be updated.  If contents are NULL, a new
 *                         entry will be added to the segment/particle
 *                         intersection list and the address of the new
 *                         structure returned to the caller in <entry>
 *                         
 *          tag1           Node tag of one of the segment endpoints
 *          tag2           Node tag of the other segment endpoint
 *          inclusionIndex Index in the local inclusion list of the 
 *                         particle intersecting the specified segment.
 *          inclusionID    ID of the particle intersecting the specified
 *                         segment.
 *
 *-------------------------------------------------------------------------*/
void SegPartListUpdate(Home_t *home, SegPartIntersect_t **entry, Tag_t *tag1,
                       Tag_t *tag2, int inclusionIndex, int inclusionID)
{
        int                nextIndex, numEntries, allocedEntries;
        SegPartIntersect_t *intersectInfo, *tmpList;

/*
 *      For threaded processes, we have to synchronize updates to the
 *      seg/particle intersection list via locks.
 */
        LOCK(&home->segpartIntersectListLock);

/*
 *      if we don't have an entry to which to add the intersect info,
 *      insert a new entry into the list.
 */
        if (*entry == (SegPartIntersect_t *)NULL) {

            numEntries = home->segpartIntersectCnt;
            allocedEntries = home->segpartIntersectListSize;
            tmpList = home->segpartIntersectList;

/*
 *          Make sure there is enough space on the list for a segment
 */
            if (numEntries >= allocedEntries) {
                allocedEntries += 100;
                tmpList = (SegPartIntersect_t *)realloc(tmpList,
                          allocedEntries * sizeof(SegPartIntersect_t));
            }

/*
 *          Get a pointer to the next free entry in the list and
 *          initialize it.  IMPORTANT: Make sure that the node tags
 *          identifying the segment are stored such that the first tag
 *          for a segment is "less than" the second as determined
 *          by OrderTags().  This is so we're consistent in the manner
 *          we sort the list and search for entries on the list later.
 */
            intersectInfo = &tmpList[numEntries];
            memset(intersectInfo, 0, sizeof(SegPartIntersect_t));

            if (OrderTags(tag1, tag2) < 0) {
                intersectInfo->tag1.domainID = tag1->domainID;
                intersectInfo->tag1.index    = tag1->index;
                intersectInfo->tag2.domainID = tag2->domainID;
                intersectInfo->tag2.index    = tag2->index;
            } else {
                intersectInfo->tag1.domainID = tag2->domainID;
                intersectInfo->tag1.index    = tag2->index;
                intersectInfo->tag2.domainID = tag1->domainID;
                intersectInfo->tag2.index    = tag1->index;
            }

            home->segpartIntersectCnt = numEntries + 1;
            home->segpartIntersectListSize = allocedEntries;
            home->segpartIntersectList = tmpList;

            *entry = intersectInfo;
        }

/*
 *      We have a pointer to the segment/particle intersection info,
 *      so just update it now
 */
        nextIndex = (*entry)->numIntersections;

        if (nextIndex < MAX_SEGPART_INTERSECTIONS)
	  {
            (*entry)->inclusionIndex[nextIndex] = inclusionIndex;
            (*entry)->inclusionID[nextIndex] = inclusionID;
            (*entry)->numIntersections += 1;
	  }
	else 
	  {
            for (int i = 0; i < nextIndex; i++)
	      {
		printf("%d %d\n",inclusionID,inclusionIndex);
	    	    
		for (int j = 0; j < (*entry)->numIntersections; j++) 
		  printf("%d %d %d\n",(*entry)->numIntersections,
			 (*entry)->inclusionIndex[j],(*entry)->inclusionID[j]);
		
		Node_t *nodetmp1;
		nodetmp1 = GetNodeFromTag(home, *tag1);
		PrintNode(nodetmp1);
		Node_t *nodetmp2;
		nodetmp2 = GetNodeFromTag(home, *tag2);
		PrintNode(nodetmp2);
	      }

	    Fatal("Segment (%d,%d)--(%d,%d) intersects"
                   " more than %d inclusions\n", tag1->domainID, tag1->index,
                   tag2->domainID, tag2->index, MAX_SEGPART_INTERSECTIONS);
	  }

        UNLOCK(&home->segpartIntersectListLock);

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     SegPartListClear()
 *      Description:  Release any memory associated with the segment/
 *                    particle intersection list and reinitialize
 *                    the associated variable
 *
 *-------------------------------------------------------------------------*/
void SegPartListClear(Home_t *home)
{
        if (home->segpartIntersectListSize > 0) {
            free(home->segpartIntersectList);
        }

        home->segpartIntersectListSize = 0;
        home->segpartIntersectCnt = 0;
        home->segpartIntersectList = (SegPartIntersect_t *)NULL;

        return;
}



/*---------------------------------------------------------------------------
 *
 *      Function:     CollisionDislocPart()
 *      Description:  Locally remesh a dislocation segment when a particle
 *                    is in the vicinity of the dislocation
 *
 *-------------------------------------------------------------------------*/
void HandleParticleCollisions(Home_t *home)
{
   
   Node_t  *node1, *node2;
                
   int thisDomain;
   thisDomain = home->myDomain;
   
   Param_t *param;
   param      = home->param;

   real8 increaseRadius = 1;

/*
 * Start looping through native nodes looking for segments to collide with particles
 */
int loopLimit = home->newNodeKeyPtr;
  // for (int i = 0; i < home->newNodeKeyPtr; i++) 
for (int i = 0; i < loopLimit; i++)
   {
      if ((node1 = home->nodeKeys[i]) == (Node_t *)NULL) continue;
      
/*
 *    Loop with all arms of node 1
 *    Only if node 1 is a discretization node
 */
      int nbrs = node1->numNbrs;

/*
 *    Add some restrictions on the segments we'll consider.
 *
 *    Discretization nodes that are neither pinned nor surface nodes.
 */
      if ((nbrs != 2) || 
          HAS_ANY_OF_CONSTRAINTS(node1->constraint,
                                 PINNED_NODE | SURFACE_NODE | PRECIPITATE_NODE)) {
         continue;
      }

      
      if ((node1->nbrTag[0].domainID != home->myDomain) ||
          (node1->nbrTag[1].domainID != home->myDomain)) {
         continue;
      }
      
      real8 pos1[3];
      pos1[X] = node1->x;
      pos1[Y] = node1->y;
      pos1[Z] = node1->z;

/*
 *    Start looping over particles now
 */
      for (int m = 0; m < home->totInclusionCount; m++) 
      {
         EInclusion_t *inclusion = &home->eshelbyInclusions[m];

/*  
 *       localMaxSeg value has to be greater than 2*minSeg+eps 
 *       so that the new created segment is not remeshed and
 *       slightly smaller than the diameter of the particle.
 */
         
	 real8 localMaxSeg = 2 * param->minSeg * 1.01;
/*
 *       Loop over the node1's arms and determine how many arms to split.
 */
         int nbrSplit = 0, armkeep[2], inckeep[2];
         armkeep[0] = -1; armkeep[1] = -1;
         inckeep[0] = -1; inckeep[1] = -1;

	 
         for (int arm12 = 0; arm12 < nbrs; arm12++)
         {
	   Node_t *node2;
	   node2 = GetNodeFromTag(home, node1->nbrTag[arm12]);
	   if (node2 == (Node_t *)NULL) continue;
	   if (!node2) continue;
	   if (!DomainOwnsSeg(home, OPCLASS_COLLISION, thisDomain, &node2->myTag)) continue;
	   
	   real8 pos2[3];
	   pos2[X] = node2->x;
	   pos2[Y] = node2->y;
	   pos2[Z] = node2->z;
           
	   inclusion->radius[X] *= increaseRadius;
	   inclusion->radius[Y] *= increaseRadius;
	   inclusion->radius[Z] *= increaseRadius;
	   
	   int ninc = 0;
	   real8 newPos1[3], newPos2[3];
	   real8 ratio[2];                  
	   int isIntersecting = IntersectionInformation(home, pos1, pos2, inclusion, 
							newPos1, newPos2, &ninc, ratio);
	   
	   inclusion->radius[X] /= increaseRadius;
	   inclusion->radius[Y] /= increaseRadius;
	   inclusion->radius[Z] /= increaseRadius;
	   
	   real8 segLen = 
	     sqrt( (newPos2[0]-newPos1[0]) * (newPos2[0]-newPos1[0]) +
		   (newPos2[1]-newPos1[1]) * (newPos2[1]-newPos1[1]) +
		   (newPos2[2]-newPos1[2]) * (newPos2[2]-newPos1[2]) );	  
		  // printf("DEBUG_COLLISION: Node(%d) SegLen=%.2e LocalMax=%.2e Intersect=%d\n",node1->myTag.index, segLen, localMaxSeg, isIntersecting); 
	   
	   if  ( (isIntersecting==1) && (segLen > localMaxSeg))
	     {
                armkeep[nbrSplit] = arm12;
                inckeep[nbrSplit] = m;
                nbrSplit += 1;

	     }
         } // end loop over arm12
                 

	 for (int isplit=0; isplit<nbrSplit; isplit++) 
	   {
/* 
 *           Split segment  
 */
	     Node_t *node2;
	     node2 = GetNodeFromTag(home, node1->nbrTag[armkeep[isplit]]);
	   
	     real8 dx, dy, dz;
	     dx = node2->x - node1->x;
	     dy = node2->y - node1->y;
	     dz = node2->z - node1->z;

	     ZImage(param, &dx, &dy, &dz);
	   
	     real8 pos2[3];
	     pos2[X] = node1->x + dx;
	     pos2[Y] = node1->y + dy;
	     pos2[Z] = node1->z + dz;
	   
	     real8 newpos[3];
	     newpos[X] = (pos2[X]+pos1[X])*0.5;
	     newpos[Y] = (pos2[Y]+pos1[Y])*0.5;
	     newpos[Z] = (pos2[Z]+pos1[Z])*0.5;
	     

// ------------------------------------------------------------
// Snap split point to (glide plane) ∩ (ellipsoid surface)
// Use the robust helper SnapPointToEllipsoidOnPlane() already defined above.
// ------------------------------------------------------------
{
    // 取当前分裂对应的析出相（按 inckeep[isplit] 最稳妥）
    EInclusion_t *incLocal = &home->eshelbyInclusions[inckeep[isplit]];

    // 取该段的滑移面：通过 node1 的该 arm 的法向给出
    int arm = armkeep[isplit];

    real8 planePt[3] = { pos1[X], pos1[Y], pos1[Z] };
    real8 planeN [3] = { node1->nx[arm], node1->ny[arm], node1->nz[arm] };

    // 先把 pin（midpoint）投影到滑移面，作为更稳定的初值
    real8 nmag = sqrt(planeN[X]*planeN[X] + planeN[Y]*planeN[Y] + planeN[Z]*planeN[Z]);
    if (nmag > 1.0e-20) {
        real8 nunit[3] = { planeN[X]/nmag, planeN[Y]/nmag, planeN[Z]/nmag };
        real8 dp = (newpos[X]-planePt[X])*nunit[X] + (newpos[Y]-planePt[Y])*nunit[Y] + (newpos[Z]-planePt[Z])*nunit[Z];
        newpos[X] -= dp*nunit[X];
        newpos[Y] -= dp*nunit[Y];
        newpos[Z] -= dp*nunit[Z];
    }

    // 计算(滑移平面 ∩ 椭球面)上最近点
    real8 snapped[3];
    if (SnapPointToEllipsoidOnPlane(home, incLocal, planePt, planeN, newpos, snapped)) {
        newpos[X] = snapped[X];
        newpos[Y] = snapped[Y];
        newpos[Z] = snapped[Z];
        FoldBox(param, &newpos[X], &newpos[Y], &newpos[Z]);
    } else {
        // fallback：至少保证还在滑移面内（不做表面贴合）
        FoldBox(param, &newpos[X], &newpos[Y], &newpos[Z]);
    }
}




	     
	     
	     real8 nodeVel[3];
	     nodeVel[X] = node2->vX;
	     nodeVel[Y] = node2->vY;
	     nodeVel[Z] = node2->vZ;

	     
	     // Split node 2
	     int splitStatus;
	     Node_t *splitNode1, *splitNode2;
	     int arm21 = GetArmID(node2, node1);

	     splitStatus = SplitNode(home,
				     OPCLASS_COLLISION,
				     node2, pos2,
				     newpos, nodeVel,
				     nodeVel, 1,
				     &arm21, 1,
				     &splitNode1,
				     &splitNode2, 0);

	     if (splitStatus != SPLIT_SUCCESS) continue;
	     // ========================================================
             // 【新增代码】 关键修正：给新节点戴上“手铐”！
             // ========================================================
             
             // SplitNode 会返回两个节点指针。
             // 通常 splitNode1 是保留的原节点，splitNode2 是新插入的节点（位于 newpos，即表面）。
             // 但为了保险，我们检查哪个节点在表面（或者直接都加上 PINNED 也不坏，
             // 因为 CheckPrecipitateCutting 会在下一时间步解开它们）。
             
             // 强制将新生成的节点标记为 PINNED
             if (splitNode1) {
                 if (splitNode1) splitNode1->constraint |= (SURFACE_NODE | PRECIPITATE_NODE);
             //    printf(">>> COLLISION TRAP: Node (%d,%d) hit particle and was PINNED!\n", 
           //splitNode1->myTag.domainID, splitNode1->myTag.index);
             }
             if (splitNode2) {
                 if (splitNode2) splitNode2->constraint |= (SURFACE_NODE | PRECIPITATE_NODE);
              //   printf(">>> COLLISION TRAP: Node (%d,%d) hit particle and was PINNED!\n", 
           //splitNode2->myTag.domainID, splitNode2->myTag.index);
             }
             
             // 只有打上了这个标签：
             // 1. 积分器(Integrator) 才会知道这个节点不能随便动。
             // 2. 你的 CheckPrecipitateCutting 才会去检查它的受力。
             // ========================================================

#ifdef FIX_PLASTIC_STRAIN
	     if (splitStatus == SPLIT_SUCCESS) {
	       real8 newPos[3];
	       
	       newPos[X] = splitNode1->x;
	       newPos[Y] = splitNode1->y;
	       newPos[Z] = splitNode1->z;
	       UpdateNodePlasticStrain(home, splitNode1, pos2, newPos);
	       
	       newPos[X] = splitNode2->x;
	       newPos[Y] = splitNode2->y;
	       newPos[Z] = splitNode2->z;
	       UpdateNodePlasticStrain(home, splitNode2, pos2, newPos);
	     }
#endif


	     MarkNodeForceObsolete(home, splitNode1);
         MarkNodeForceObsolete(home, splitNode2);
	   }	 
	 
      } // end loop over inclusions
   } // loop over all nodes
}

#endif // ESHELBY
/*---------------------------------------------------------------------------
 *
 * Function:     CheckPrecipitateCutting
 * Description:  Checks if pinned nodes have enough force to cut through
 * the particle.
 * * Logic:        1. Reads RAW force from ParaDiS (Magnitude ~10^8).
 * 2. Converts it to NEWTONS by multiplying by b^2.
 * 3. Compares with Threshold (in Newtons).
 * 4. Unpins node if Force > Threshold.
 *
 *-------------------------------------------------------------------------*/
void CheckPrecipitateCutting(Home_t *home)
{
#ifndef ESHELBY
    (void)home;
    return;
#else
    Param_t *param = home->param;

    // ---------------------------------------------------------
    // 1) Critical cutting force threshold (Newtons)
    //    (default can be overridden via environment variable
    //     PRECIPITATE_STRENGTH)
    // ---------------------------------------------------------
    static int    initialized            = 0;
    static real8  criticalForceThreshold = 1.0e-9;  // N

    if (!initialized) {
        char *envVal = getenv("PRECIPITATE_STRENGTH");
        if (envVal != (char *)NULL) {
            real8 tmp = atof(envVal);
            if (tmp > 0.0) {
                criticalForceThreshold = tmp;
            }
        }
        initialized = 1;

        printf("[Precipitate] Critical cutting force threshold = %e N\n",
               criticalForceThreshold);
    }

    // Convert internal nodal force units (Pa*b^2) to Newtons.
    // Geometry length unit is b; burgMag provides b in meters:
    // F[N] = F[Pa*b^2] * (burgMag[m])^2
    real8 unit_conversion = param->burgMag * param->burgMag;

    // ---------------------------------------------------------
    // 2) Surface-contact enforcement parameters (dimensionless)
    //    phi = (x'/a)^2 + (y'/b)^2 + (z'/c)^2
    // ---------------------------------------------------------
    const real8 gapStick   = 0.02;  // if |phi-1| <= gapStick => enforce v·n = 0
    const real8 gapRelease = 0.10;  // if phi > 1+gapRelease => release contact

    for (int i = 0; i < home->newNodeKeyPtr; i++) {

        Node_t *node = home->nodeKeys[i];
        if (node == (Node_t *)NULL) {
            continue;
        }

        if (!HAS_ANY_OF_CONSTRAINTS(node->constraint, PRECIPITATE_NODE)) {
            continue;
        }

        // Keep these nodes in the partial force/velocity update set
        MarkNodeForceObsolete(home, node);

        if (home->totInclusionCount <= 0) {
            node->constraint &= ~(PRECIPITATE_NODE | SURFACE_NODE);
            continue;
        }

        // -----------------------------------------------------
        // 3) Find nearest inclusion (min |phi-1|)
        // -----------------------------------------------------
        int   bestInc = -1;
        real8 bestPhi = 1.0e30;
        real8 bestGap = 1.0e30;

        real8 bestU1[3], bestU2[3], bestU3[3];
        real8 bestLoc[3];
        real8 bestRad[3];

        for (int m = 0; m < home->totInclusionCount; m++) {

            EInclusion_t *inc = &home->eshelbyInclusions[m];

            real8 ctr[3];
            ctr[X] = inc->position[X];
            ctr[Y] = inc->position[Y];
            ctr[Z] = inc->position[Z];

            PBCPOSITION(param, node->x, node->y, node->z,
                        &ctr[X], &ctr[Y], &ctr[Z]);

            real8 u1[3], u2[3], u3[3];
            u1[X] = inc->rotation[0][0];
            u1[Y] = inc->rotation[0][1];
            u1[Z] = inc->rotation[0][2];
            u2[X] = inc->rotation[1][0];
            u2[Y] = inc->rotation[1][1];
            u2[Z] = inc->rotation[1][2];
            u3[X] = u1[Y]*u2[Z] - u1[Z]*u2[Y];
            u3[Y] = u1[Z]*u2[X] - u1[X]*u2[Z];
            u3[Z] = u1[X]*u2[Y] - u1[Y]*u2[X];

            real8 d[3];
            d[X] = node->x - ctr[X];
            d[Y] = node->y - ctr[Y];
            d[Z] = node->z - ctr[Z];

            real8 loc[3];
            loc[X] = d[X]*u1[X] + d[Y]*u1[Y] + d[Z]*u1[Z];
            loc[Y] = d[X]*u2[X] + d[Y]*u2[Y] + d[Z]*u2[Z];
            loc[Z] = d[X]*u3[X] + d[Y]*u3[Y] + d[Z]*u3[Z];

            real8 aR = inc->radius[X];
            real8 bR = inc->radius[Y];
            real8 cR = inc->radius[Z];

            if ((aR <= 0.0) || (bR <= 0.0) || (cR <= 0.0)) {
                continue;
            }

            real8 phi = (loc[X]*loc[X])/(aR*aR) +
                        (loc[Y]*loc[Y])/(bR*bR) +
                        (loc[Z]*loc[Z])/(cR*cR);

            real8 gap = fabs(phi - 1.0);

            if (gap < bestGap) {
                bestGap = gap;
                bestPhi = phi;
                bestInc = m;

                bestRad[X] = aR;
                bestRad[Y] = bR;
                bestRad[Z] = cR;

                bestLoc[X] = loc[X];
                bestLoc[Y] = loc[Y];
                bestLoc[Z] = loc[Z];

                bestU1[X] = u1[X]; bestU1[Y] = u1[Y]; bestU1[Z] = u1[Z];
                bestU2[X] = u2[X]; bestU2[Y] = u2[Y]; bestU2[Z] = u2[Z];
                bestU3[X] = u3[X]; bestU3[Y] = u3[Y]; bestU3[Z] = u3[Z];
            }
        }

        // If the node drifted away or we cannot associate it with any inclusion
        if ((bestInc < 0) || (bestGap > 0.5)) {
            node->constraint &= ~(PRECIPITATE_NODE | SURFACE_NODE);
            continue;
        }

        // -----------------------------------------------------
        // Compute outward normal from grad(phi) in global coords
        // -----------------------------------------------------
        real8 gradLoc[3];
        gradLoc[X] = bestLoc[X] / (bestRad[X]*bestRad[X]);
        gradLoc[Y] = bestLoc[Y] / (bestRad[Y]*bestRad[Y]);
        gradLoc[Z] = bestLoc[Z] / (bestRad[Z]*bestRad[Z]);

        real8 n[3];
        n[X] = gradLoc[X]*bestU1[X] + gradLoc[Y]*bestU2[X] + gradLoc[Z]*bestU3[X];
        n[Y] = gradLoc[X]*bestU1[Y] + gradLoc[Y]*bestU2[Y] + gradLoc[Z]*bestU3[Y];
        n[Z] = gradLoc[X]*bestU1[Z] + gradLoc[Y]*bestU2[Z] + gradLoc[Z]*bestU3[Z];

        real8 nMag = sqrt(n[X]*n[X] + n[Y]*n[Y] + n[Z]*n[Z]);
        if (nMag > 1.0e-20) {
            n[X] /= nMag;
            n[Y] /= nMag;
            n[Z] /= nMag;
        }

        // -----------------------------------------------------
// Surface stage with NO cross-slip:
// enforce v lies in (glide plane) ∩ (surface tangent plane)
// i.e., v || (g x n), where g = glide plane normal, n = surface normal
// -----------------------------------------------------
{
    // 1) Build an effective glide normal g from node arms
    real8 g[3] = {0.0, 0.0, 0.0};
    int gCount = 0;

    for (int a = 0; a < node->numNbrs; a++) {
        real8 ng[3] = { node->nx[a], node->ny[a], node->nz[a] };
        real8 ngm = sqrt(ng[X]*ng[X] + ng[Y]*ng[Y] + ng[Z]*ng[Z]);
        if (ngm < 1.0e-20) continue;
        ng[X] /= ngm; ng[Y] /= ngm; ng[Z] /= ngm;

        // align sign to avoid cancellation
        if (gCount > 0) {
            real8 dot = ng[X]*g[X] + ng[Y]*g[Y] + ng[Z]*g[Z];
            if (dot < 0.0) { ng[X] = -ng[X]; ng[Y] = -ng[Y]; ng[Z] = -ng[Z]; }
        }
        g[X] += ng[X]; g[Y] += ng[Y]; g[Z] += ng[Z];
        gCount++;
    }

    real8 gm = sqrt(g[X]*g[X] + g[Y]*g[Y] + g[Z]*g[Z]);
    if (gm > 1.0e-20) {
        g[X] /= gm; g[Y] /= gm; g[Z] /= gm;

        // 2) Intersection direction t = g x n
        real8 t[3];
        t[X] = g[Y]*n[Z] - g[Z]*n[Y];
        t[Y] = g[Z]*n[X] - g[X]*n[Z];
        t[Z] = g[X]*n[Y] - g[Y]*n[X];

        real8 tm2 = t[X]*t[X] + t[Y]*t[Y] + t[Z]*t[Z];

        if (tm2 > 1.0e-30) {
            real8 invtm = 1.0 / sqrt(tm2);
            t[X] *= invtm; t[Y] *= invtm; t[Z] *= invtm;

            // project current v onto t
            real8 vt = node->vX*t[X] + node->vY*t[Y] + node->vZ*t[Z];
            node->vX = vt*t[X];
            node->vY = vt*t[Y];
            node->vZ = vt*t[Z];

        } else {
            // g parallel to n: just keep glide (no out-of-plane)
            real8 vg = node->vX*g[X] + node->vY*g[Y] + node->vZ*g[Z];
            node->vX -= vg*g[X];
            node->vY -= vg*g[Y];
            node->vZ -= vg*g[Z];

            // and prevent inward penetration if still inside
            real8 vn2 = node->vX*n[X] + node->vY*n[Y] + node->vZ*n[Z];
            if ((bestPhi < 1.0) && (vn2 < 0.0)) {
                node->vX -= vn2*n[X];
                node->vY -= vn2*n[Y];
                node->vZ -= vn2*n[Z];
            }
        }
    } else {
        // If we cannot determine glide normal, fall back to surface tangential only
        real8 vn = node->vX*n[X] + node->vY*n[Y] + node->vZ*n[Z];
        if (bestGap <= gapStick) {
            node->vX -= vn*n[X];
            node->vY -= vn*n[Y];
            node->vZ -= vn*n[Z];
        } else if ((bestPhi < 1.0) && (vn < 0.0)) {
            node->vX -= vn*n[X];
            node->vY -= vn*n[Y];
            node->vZ -= vn*n[Z];
        }
    }
}

        // Release if separated sufficiently (bypass complete)
        if (bestPhi > 1.0 + gapRelease) {
            node->constraint &= ~(PRECIPITATE_NODE | SURFACE_NODE);
            continue;
        }

        // -----------------------------------------------------
        // Cutting criterion (based on nodal force magnitude)
        // -----------------------------------------------------
        real8 fx = node->fX;
        real8 fy = node->fY;
        real8 fz = node->fZ;

        real8 forceMagInternal = sqrt(fx*fx + fy*fy + fz*fz);
        real8 physicalForce    = forceMagInternal * unit_conversion;

        if (physicalForce > criticalForceThreshold) {
            node->constraint &= ~(PRECIPITATE_NODE | SURFACE_NODE);

            printf(">>> CUTTING EVENT: Node (%d,%d) released from precipitate contact.\n"
       "    |f_internal| = %e , b^2 = %e , Force = %e N, Threshold = %e N\n",
       node->myTag.domainID, node->myTag.index,
       forceMagInternal, unit_conversion, physicalForce, criticalForceThreshold);
        }
    }
#endif
}