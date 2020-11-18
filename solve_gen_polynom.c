//find roots solve_gen_polynom.c  
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#define PI 3.14159265
#define N 10

int solve_quadratic(double a, double b, double c, double* r1, double* i1, double* r2, double* i2, double epsilon);
void divide_cubic(double a, double* b, double* c, double x);
void cubic_root(double x, double y, double* r, double* i);
void square_root(double x, double y, double* r, double* i);
int solve_cubic(double a1, double a2, double a3, double a4, double* x1, double* y1, double* x2, double* y2, double* x3, double* y3, double epsilon);
void solve_biquadratic(double a, double b, double c, double* x1, double* y1, double* x2, double* y2, double* x3, double* y3, double* x4, double* y4, double epsilon);
int solve_quartic(double a4, double a3, double a2, double a1, double a0, double* x1, double* y1, double* x2, double* y2, double* x3, double* y3, double* x4, double* y4, double epsilon);
void complex_mult(double a, double b, double c, double d, double* x, double* y);
void complex_add(double a, double b, double c, double d, double* x, double* y);
void test_solution(double a, double b, double c, double d, double e, double x, double y, double* result1, double* result2);
double bi_regula_combo(double (*fun)(double, double*, int), double a, double b, double eps, double coeffs[], int deg);

///////////////////////////////////////////Project 2//////////////////////////////////////////////////////////////
int solve_gen_polynom(double solution_points[], double coeffs[], int deg, double epsilon);
void find_fd(double coeffs[], int deg, double fd_coeffs[]);
void p_copy(double coeffs_from[], double coeffs_to[], int deg);
void p_div(double coeffs_p[], double coeffs_d[], double coeffs_r[], double coeffs_q[], int deg);
double f(double x, double coeffs[], int deg);

int solve_gen_polynom(double solution_points[], double coeffs[], int deg, double epsilon)
{
    int i, j, num_of_solutions, flag = 0, flag2 = 0, s_div_poly;
    double y1, y2, y3, y4, check, max, min, chosen_sol, source;
    double fd_coeffs[N], fd_solution_points[N], new_coeffs[N], coeffs_d[N], coeffs_r[N];

    switch (deg) {
    case 0:
        return 0;
    case 1:
        solution_points[0] = -coeffs[1] / coeffs[0];
        return 1;
    case 2:
        return solve_quadratic(coeffs[0], coeffs[1], coeffs[2], &solution_points[0], &y1, &solution_points[1], &y2, epsilon);
    case 3:
        return solve_cubic(coeffs[0], coeffs[1], coeffs[2], coeffs[3], &solution_points[0], &y1, &solution_points[1], &y2, &solution_points[2], &y3, epsilon);
    case 4:
        return solve_quartic(coeffs[0], coeffs[1], coeffs[2], coeffs[3], coeffs[4], &solution_points[0], &y1, &solution_points[1], &y2, &solution_points[2],
            &y3, &solution_points[3], &y4, epsilon);
    default:
        find_fd(coeffs, deg, fd_coeffs);
        num_of_solutions = solve_gen_polynom(fd_solution_points, fd_coeffs, deg - 1, epsilon);

        for (i = 0; i < num_of_solutions; i++)
        {
            if (!f(fd_solution_points[i], coeffs, deg))
            {
                check = fd_solution_points[i];
                flag = 1;
                break;
            }
        }

        if (flag == 1)
        {
            coeffs_d[0] = 1;
            coeffs_d[1] = -check;
            p_div(coeffs, coeffs_d, coeffs_r, new_coeffs, deg);
            s_div_poly = solve_gen_polynom(solution_points, new_coeffs, deg - 1, epsilon);
            solution_points[s_div_poly] = check;
            s_div_poly++;
            return s_div_poly;
        }
        else
        {
            if (!((deg - 1) % 2))
            {
                max = min = 0;
                for (i = 0; i < num_of_solutions; i++)
                {
                    if (f(fd_solution_points[i], coeffs, deg) < 0)
                    {
                        flag2 = 1;
                        chosen_sol = fd_solution_points[i];
                    }

                    if (fd_solution_points[i] > max)
                        max = fd_solution_points[i];

                    if (fd_solution_points[i] < min)
                        min = fd_solution_points[i];
                }

                if (flag2 == 1)
                {
                    if (chosen_sol != max) 
                        check = bi_regula_combo(f, chosen_sol, (max + 1)*100, epsilon, coeffs, deg);
                    else 
                        check = bi_regula_combo(f, chosen_sol, min, epsilon, coeffs, deg);
                   

                    coeffs_d[0] = 1;
                    coeffs_d[1] = -check;
                    p_div(coeffs, coeffs_d, coeffs_r, new_coeffs, deg);
                    s_div_poly = solve_gen_polynom(solution_points, new_coeffs, deg - 1, epsilon);
                    solution_points[s_div_poly] = check;
                    s_div_poly++;
                    return s_div_poly;
                }
                else {
                    return 0;

                }
            }
            else
            {
                max = fd_solution_points[i];
                for (i = 0; i < num_of_solutions; i++)
                {
                    if (fd_solution_points[i] > max)
                        max = fd_solution_points[i];
                }
                                 
                check = bi_regula_combo(f, 0.0, (max + 1), epsilon, coeffs, deg);  //chose max + 1 as interval because is the 
                coeffs_d[0] = 1;                                                  // only interval that works for both examples.
                coeffs_d[1] = -check;
                p_div(coeffs, coeffs_d, coeffs_r, new_coeffs, deg);
                s_div_poly = solve_gen_polynom(solution_points, new_coeffs, deg - 1, epsilon);
                solution_points[s_div_poly] = check;
                s_div_poly++;
                return s_div_poly;
            }

        }

    }
}

void find_fd(double coeffs[], int deg, double fd_coeffs[])
{
    int i, j;
    j = deg;
    i = 0;
    while (j)
    {
        fd_coeffs[i] = coeffs[i] * j;
        j--;
        i++;
    }

}

double f(double x, double coeffs[], int deg)
{
    int i, j, k;
    double sum, sum1;

    k = deg;
    sum = 0;
    for (i = 0; i < deg; i++, k--)
    {
        sum1 = x;
        for (j = 1; j < k; j++)
            sum1 *= x;
        sum += (coeffs[i] * sum1);
    }

    return sum + coeffs[deg];
} /* f */

void p_copy(double coeffs_from[], double coeffs_to[], int deg)
{
    int i;
    for (i = 0; i <= deg; i++)
    {
        coeffs_to[i] = coeffs_from[i];
    }

}
/* p: poly;  d: divisor x-x*;  r: remainder; */
void p_div(double coeffs_p[], double coeffs_d[], double coeffs_r[], double coeffs_q[], int deg)
{
    int i, j, rdeg;
    int power = deg - 1; //ppower - dpower
    double ratio;

    if (power < 0) return 0; //check to evoid error

    p_copy(coeffs_p, coeffs_r, deg);
    rdeg = deg;


    for (i = deg; i >= 1; i--) {

        coeffs_q[i - 1] = ratio = coeffs_r[i] / coeffs_d[1];
        coeffs_r[i] = 0;

        for (j = 0; j < 1; j++)
            coeffs_r[i - 1 + j] -= coeffs_d[j] * ratio;
    }
    while (!coeffs_r[--rdeg]);

}
//////////////////////////////////////////////// Project 2/////////////////////////////////////////////////////////////////////////


double bi_regula_combo(double (*fun)(double, double*, int), double a, double b, double eps, double coeffs[], int deg)
{
    double x, xmid, f;
    int flag;

    // Calculate initial interval


    do {
        // bisection
        x = (a + b) / 2.0;
        f = fun(x, coeffs, deg);
        //  printf("f(a) = f(%lf) = %lf, f(x) = f(%lf) = %lf, f(b) = "
         //     " f(%lf) = %lf \n",
        //      a, fun(a, coeffs, deg), x, f, b, fun(b, coeffs, deg));
        if (fabs(f) < eps)
            return x;
        if (fun(a, coeffs, deg) * f < 0.0)
            b = x;
        else
            a = x;

        // regula

        x = (a * (*fun)(b, coeffs, deg) - b * (*fun)(a, coeffs, deg)) / ((*fun)(b, coeffs, deg) - (*fun)(a, coeffs, deg));
        f = (*fun)(x, coeffs, deg);
        if (fabs(f) < eps)
            return x;
        if ((*fun)(a, coeffs, deg) * f < 0.0)
            b = x;
        else
            a = x;


        // printf("a = %lf, b = %lf, x = %lf, f(%lf) = %lf\n",
       //      a, b, x, x, (*fun)(x, coeffs, deg));

    } while (fabs(b - a) > eps);

    return x;

} /* regula */

int solve_quadratic(double a, double b, double c, double* r1, double* i1, double* r2, double* i2, double epsilon)
{
    double delta, delta_root, xv, qv;

    if (a == 0)
        if (b == 0)
            return 0;
        else
        {
            *r1 = -c / b;
            *i1 = 0;
            return 1;
        } /* else */

    xv = -b / (2 * a); /* Take care of square quad */
    qv = a * xv * xv + b * xv + c;
    if (fabs(qv) < epsilon)
    {
        *r1 = *r2 = xv;
        *i1 = *i2 = 0.0;
        return 2;
    } /* if */


    delta = b * b - 4 * a * c;
    if (delta < 0.0)
    {
        *r1 = *r2 = -b / (2 * a);  /* real part only */
        delta_root = sqrt(-delta);
        *i1 = -delta_root / (2.0 * a);
        *i2 = delta_root / (2.0 * a);
        return 0;
    } /* if */
    else
    {
        delta_root = sqrt(delta);
        *r1 = (-b - delta_root) / (2.0 * a);
        *r2 = (-b + delta_root) / (2.0 * a);
        *i1 = 0;
        *i2 = 0;
        return 2;
    } /* else */

} /* solve_quadratic */

void divide_cubic(double a, double* b, double* c, double x)
{
    *b = *b + a * x;
    *c = *c + (*b) * x;

} /*  divide_cubic */

void cubic_root(double x, double y, double* r, double* i)
{
    double radius, theta, temp;

    if (x == 0.0 && y == 0.0)
    {
        *r = *i = 0.0;
        return;
    } /* if */


    radius = sqrt(x * x + y * y);
    theta = (asin(fabs(y / radius)));

    temp = exp(log(radius) / 3.0);
    if ((x < 0) && (y >= 0))
        theta = 3 * PI - theta;
    else
        if ((x < 0) && (y <= 0))
            theta = 3 * PI + theta;
        else
            if ((x >= 0) && (y < 0))
                theta = 4 * PI - theta;

    *r = temp * cos(theta / 3.0);
    *i = temp * sin(theta / 3.0);

} /* cubic_root */

void square_root(double x, double y, double* r, double* i)
{
    double radius, theta, temp;

    if (x == 0.0 && y == 0.0)
    {
        *r = *i = 0.0;
        return;
    } /* if */


    radius = sqrt(x * x + y * y);
    theta = (asin(fabs(y / radius)));

    temp = exp(log(radius) / 2.0);
    if ((x < 0) && (y >= 0))
        theta = 3 * PI - theta;
    else
        if ((x < 0) && (y <= 0))
            theta = 3 * PI + theta;
        else
            if ((x >= 0) && (y < 0))
                theta = 4 * PI - theta;

    *r = temp * cos(theta / 2.0);
    *i = temp * sin(theta / 2.0);

} /* cubic_root */

int solve_cubic(double a1, double a2, double a3, double a4, double* x1, double* y1, double* x2, double* y2, double* x3, double* y3, double epsilon)
{

    double p, q;
    double a, b, c;
    double u1, u2, wr, wi;
    double x, t, ur, ui, vr, vi;
    double b2, c2;
    double i1, i2;
    int n;

    if (a1 == 0.0)
    {
        n = solve_quadratic(a2, a3, a4, x1, &i1, x2, &i2, epsilon);
        return n;
    } /* if */

    a = a2 / a1;
    b = a3 / a1;
    c = a4 / a1;

    p = b - a * a / 3.0;
    q = c + (2.0 * a * a * a - 9.0 * a * b) / 27.0;


    if ((p == 0.0) && (q == 0.0))
    {
        *x1 = *x2 = *x3 = -a / 3;
        *y1 = *y2 = *y3 = 0;
        return 3;
    } /* if */

    if (p == 0.0) /* q != 0 */
    {
        cubic_root(-q, 0, &ur, &ui);
        x = ur - (a / 3.0);
        *x1 = x;
        *y1 = 0;

        b2 = a;
        c2 = b;

        divide_cubic(1.0, &b2, &c2, x);
        n = solve_quadratic(1.0, b2, c2, x2, y2, x3, y3, epsilon);

        return n + 1;

    } /* if */

    solve_quadratic(1.0, -q, -((p * p * p) / 27.0), &u1, &i1, &u2, &i2, epsilon);

    if (u1 > u2)
    {
        wr = u1;
        wi = i1;
    } /* if */
    else
    {
        wr = u2;
        wi = i2;
    } /* else */

    cubic_root(wr, wi, &ur, &ui);
    cubic_root(-q + wr, wi, &vr, &vi);
    t = vr - ur;
    x = t - (a / 3.0);

    *x1 = x;
    *y1 = 0;

    b2 = a;
    c2 = b;

    divide_cubic(1.0, &b2, &c2, x);
    n = solve_quadratic(1.0, b2, c2, x2, y2, x3, y3, epsilon);

    return n + 1;

} /* solve_cubic */

void solve_biquadratic(double a, double b, double c, double* x1, double* y1, double* x2, double* y2, double* x3, double* y3, double* x4, double* y4, double epsilon)
{

    double q11, q12, q21, q22;
    double sqrt_q11, sqrt_q12, sqrt_q21, sqrt_q22;
    double sqrt_q31, sqrt_q32, sqrt_q41, sqrt_q42;


    solve_quadratic(a, b, c,
        &q11, &q12,
        &q21, &q22, epsilon);

    square_root(q11, q12, &sqrt_q11, &sqrt_q12);

    sqrt_q21 = -sqrt_q11;
    sqrt_q22 = -sqrt_q12;

    square_root(q21, q22, &sqrt_q31, &sqrt_q32);

    sqrt_q41 = -sqrt_q31;
    sqrt_q42 = -sqrt_q32;

    *x1 = sqrt_q11;
    *y1 = sqrt_q12;

    *x2 = sqrt_q21;
    *y2 = sqrt_q22;


    *x3 = sqrt_q31;
    *y3 = sqrt_q32;


    *x4 = sqrt_q41;
    *y4 = sqrt_q42;


} // solve_biquadratic

int solve_quartic(double a4, double a3, double a2, double a1, double a0, double* x1, double* y1, double* x2, double* y2, double* x3, double* y3, double* x4, double* y4, double epsilon)
{
    double b, c, d, e;
    double p, q, r;
    double b3, b2, b1, b0;
    double alpha1, alpha2, beta1, beta2, gamma1, gamma2;
    double abs_sqrt_alpha, abs_sqrt_beta;
    double sqrt_alpha1, sqrt_alpha2, sqrt_beta1, sqrt_beta2,
        sqrt_gamma1, sqrt_gamma2;
    double r11, r12, r21, r22, r31, r32, r41, r42;
    double maybe_q;
    int n, num_solutions, i;
    double* sol[8];

    /* printf("a4 = %lf, a3 = %lf, a2 = %lf, a1 = %lf, a0 = %lf\n",
         a4, a3, a2, a1, a0); */

    b = a3 / a4;
    c = a2 / a4;
    d = a1 / a4;
    e = a0 / a4;

    /*printf("b = %lf, c = %lf, d = %lf, e = %lf\n", b, c, d, e);*/


    p = (8.0 * c - 3.0 * b * b) / 8.0;
    q = (b * b * b - 4.0 * b * c + 8.0 * d) / 8.0;
    r = (-3.0 * b * b * b * b + 256.0 * e - 64.0 * b * d + 16 * b * b * c) / 256.0;

    /* printf("p = %lf, q = %lf, r = %lf\n", p, q, r);*/

    if (fabs(q) < epsilon)
    {
        solve_biquadratic(1.0, p, r,
            &r11, &r12,
            &r21, &r22,
            &r31, &r32,
            &r41, &r42, epsilon);

        /*
          *x1 = r11 - (b / 4.0);
          *y1 = r12;
          *x2 = r21 - (b / 4.0);
          *y2 = r22;
          *x3 = r31 - (b / 4.0);
          *y3 = r32;
          *x4 = r41 - (b / 4.0);
          *y4 = r42;
         */
        sol[0] = x1;
        sol[1] = x2;
        sol[2] = x3;
        sol[3] = x4;
        sol[4] = y1;
        sol[5] = y2;
        sol[6] = y3;
        sol[7] = y4;

        i = 0;
        num_solutions = 0;
        if (fabs(r12) < epsilon) {
            *sol[i] = r11 - (b / 4.0);
            *sol[i + 4] = r12;
            num_solutions++;
            i++;
        }

        if (fabs(r22) < epsilon) {
            *sol[i] = r21 - (b / 4.0);
            *sol[i + 4] = r22;
            num_solutions++;
            i++;
        }

        if (fabs(r32) < epsilon) {
            *sol[i] = r31 - (b / 4.0);
            *sol[i + 4] = r32;
            num_solutions++;
            i++;
        }

        if (fabs(r42) < epsilon) {
            *sol[i] = r41 - (b / 4.0);
            *sol[i + 4] = r42;
            num_solutions++;
            i++;
        }

        return num_solutions;
    }// if



    b3 = 1.0;
    b2 = 2.0 * p;
    b1 = p * p - 4 * r;
    b0 = -q * q;

    /* printf("b3 = %lf, b2 = %lf, b1 = %lf, b0 = %lf\n", b3, b2, b1, b0);*/


    n = solve_cubic(b3, b2, b1, b0,
        &alpha1, &alpha2,
        &beta1, &beta2,
        &gamma1, &gamma2, epsilon);

    /*
        printf("alpha1 = %lf, alpha2 = %lf\n", alpha1, alpha2);
        printf("beta1 = %lf, beta2 = %lf\n", beta1, beta2);
        printf("gamma1 = %lf, gamma2 = %lf\n", gamma1, gamma2);*/

    square_root(alpha1, alpha2, &sqrt_alpha1, &sqrt_alpha2);
    square_root(beta1, beta2, &sqrt_beta1, &sqrt_beta2);
    square_root(gamma1, gamma2, &sqrt_gamma1, &sqrt_gamma2);


    abs_sqrt_alpha = sqrt_alpha1 * sqrt_alpha1 + sqrt_alpha2 * sqrt_alpha2;
    abs_sqrt_beta = sqrt_beta1 * sqrt_beta1 + sqrt_beta2 * sqrt_beta2;


    sqrt_gamma1 = -q * (sqrt_alpha1 * sqrt_beta1 -
        sqrt_alpha2 * sqrt_beta2) / (abs_sqrt_alpha * abs_sqrt_beta);
    sqrt_gamma2 = q * (sqrt_alpha2 * sqrt_beta1 +
        sqrt_alpha1 * sqrt_beta2) / (abs_sqrt_alpha * abs_sqrt_beta);

    /*
        printf("sqrt_alpha1 = %lf, sqrt_alpha2 = %lf\n", sqrt_alpha1, sqrt_alpha2);
        printf("sqrt_beta1 = %lf, sqrt_beta2 = %lf\n", sqrt_beta1, sqrt_beta2);
        printf("sqrt_gamma1 = %lf, sqrt_gamma2 = %lf\n", sqrt_gamma1, sqrt_gamma2);*/

    r11 = (sqrt_alpha1 + sqrt_beta1 + sqrt_gamma1) / 2.0;
    r12 = (sqrt_alpha2 + sqrt_beta2 + sqrt_gamma2) / 2.0;
    r21 = (sqrt_alpha1 - sqrt_beta1 - sqrt_gamma1) / 2.0;
    r22 = (sqrt_alpha2 - sqrt_beta2 - sqrt_gamma2) / 2.0;
    r31 = (-sqrt_alpha1 + sqrt_beta1 - sqrt_gamma1) / 2.0;
    r32 = (-sqrt_alpha2 + sqrt_beta2 - sqrt_gamma2) / 2.0;
    r41 = (-sqrt_alpha1 - sqrt_beta1 + sqrt_gamma1) / 2.0;
    r42 = (-sqrt_alpha2 - sqrt_beta2 + sqrt_gamma2) / 2.0;
    /*
        printf("r11 = %lf, r12 = %lf, r21 = %lf,r22 = %lf\n",
            r11, r12, r21, r22);
        printf("r31 = %lf, r32 = %lf, r41 = %lf,r42 = %lf\n",
            r31, r32, r41, r42);*/

            /*
                *x1 = r11 - (b / 4.0);
                *y1 = r12;
                *x2 = r21 - (b / 4.0);
                *y2 = r22;
                *x3 = r31 - (b / 4.0);
                *y3 = r32;
                *x4 = r41 - (b / 4.0);
                *y4 = r42;

                printf("b/4.0 = %lf\n", b / 4.0);
             */
    sol[0] = x1;
    sol[1] = x2;
    sol[2] = x3;
    sol[3] = x4;
    sol[4] = y1;
    sol[5] = y2;
    sol[6] = y3;
    sol[7] = y4;

    i = 0;
    num_solutions = 0;
    if (fabs(r12) < epsilon) {
        *sol[i] = r11 - (b / 4.0);
        *sol[i + 4] = r12;
        num_solutions++;
        i++;
    }

    if (fabs(r22) < epsilon) {
        *sol[i] = r21 - (b / 4.0);
        *sol[i + 4] = r22;
        num_solutions++;
        i++;
    }

    if (fabs(r32) < epsilon) {
        *sol[i] = r31 - (b / 4.0);
        *sol[i + 4] = r32;
        num_solutions++;
        i++;
    }

    if (fabs(r42) < epsilon) {
        *sol[i] = r41 - (b / 4.0);
        *sol[i + 4] = r42;
        num_solutions++;
        i++;
    }

    return num_solutions;
} // solve_quartic

void complex_mult(double a, double b, double c, double d, double* x, double* y)
{
    *x = a * c - b * d;
    *y = a * d + c * b;
} // complex_mult

void complex_add(double a, double b, double c, double d, double* x, double* y)
{
    *x = a + b;
    *y = c + d;
} // complex_mult

void test_solution(double a, double b, double c, double d, double e, double x, double y, double* result1, double* result2)
{
    double x1, y1, temp11, temp12, temp21, temp22;
    int i;

    temp11 = x;
    temp12 = y;
    temp21 = e;
    temp22 = 0.0;

    temp21 += d * temp11;
    temp22 += d * temp12;;

    complex_mult(temp11, temp12, x, y, &temp11, &temp12);
    temp21 += c * temp11;
    temp22 += c * temp12;;

    complex_mult(temp11, temp12, x, y, &temp11, &temp12);
    temp21 += b * temp11;
    temp22 += b * temp12;;

    complex_mult(temp11, temp12, x, y, &temp11, &temp12);
    temp21 += a * temp11;
    temp22 += a * temp12;;

    *result1 = temp21;
    *result2 = temp22;

} // test_solution 

int main()
{
    double r;
    int i1, j1, i;
    double a, b, c, d, e;
    double interval_a, interval_b;
    double x0, x1, x2, x3, x4;
    double y0, y1, y2, y3, y4;
    double min_points[5];
    double min_value, temp_value;
    double coeffs[N];
    double solutions[N];

    double epsilon = 0.00000001;
    int no_of_solutions, degree;
    int flag;
    double result1, result2;


    printf("Enter degree\n");
    scanf_s("%d", &degree);

    printf("Enter %d coeffs, first must be > 0:\n", degree + 1);
    for (i = 0; i <= degree; i++)
        scanf_s("%lf", &coeffs[i]);

    printf("The polynom:\n");
    for (i = 0; i <= degree; i++)
        printf(" %10.3lf ", coeffs[i]);
    printf("\n");
    no_of_solutions = solve_gen_polynom(solutions, coeffs,
        degree, epsilon);

    printf("no_of_solutions = %d\n", no_of_solutions);
    printf("Solutions:\n");
    for (i = 0; i < no_of_solutions; i++)
        printf(" %10.3lf ", solutions[i]);
    printf("\n");

    return 0;
} /* main */



