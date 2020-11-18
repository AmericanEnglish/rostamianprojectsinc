/* The cut-off function for demo-square-with-hole.c is
 *     phi(x,y) = (1-x^2)*(1-y^2)*(x^2 + y^2 - 1/4).
 * The derivatives of phi where calculated symbolically
 * in Maple.  The result was translated to C code also in Maple.
*/
static void my_phi(struct Neural_Net_PDE *nn, double x, double y)
{
	double t1 = x * x;
	double t2 = -t1 + 0.1e1;
	double t4 = y * y;
	double t5 = -t4 + 0.1e1;
	double t6 = t1 + t4 - 0.1e1 / 0.4e1;
	double t7 = t5 * t6;
	double t9 = x * t5;
	double t11 = t2 * t5;
	double t15 = t2 * y;
	double t23 = 0.2e1 * t11;

	nn->phi[0][0] = t2 * t7;
	nn->phi[1][0] = 0.2e1 * (t11 * x - t9 * t6);
	nn->phi[0][1] = 0.2e1 * (t11 * y - t15 * t6);
	nn->phi[2][0] = -0.8e1 * t1 * t5 + t23 - 0.2e1 * t7;
	nn->phi[1][1] = 0.4e1 * (x * y * t6 - t15 * x - t9 * y);
	nn->phi[0][2] = -0.8e1 * t2 * t4 - 0.2e1 * t2 * t6 + t23;
}

/* For the exact solution we pick u(x,t) = (2-x)*phi(x,y),
 * and then we "reverse engineer" the PDE by picking
 * f = u_xx + u_yy.  The calculation of f is not difficult
 * by hand, but I did it in Maple to avoid risk of errors.
*/
static double my_pde(double x, double y, double u,
		double u_x, double u_y,
		double u_xx, double u_xy, double u_yy)
{
	double t1 = x * x;
	double t2 = t1 * t1;
	double t6 = y * y;
	double t16 = t6 * t6;
	double f = 0.2e1 * t2 * x + 0.4e1 * t2
		+ (0.64e2 * t6 - 0.49e2) * t1 * x / 0.2e1
		+ (0.96e2 * t6 - 0.66e2) * t1 / 0.2e1
		+ (0.12e2 * t16 - 0.51e2 * t6 + 0.20e2) * x / 0.2e1
		+ 0.4e1 * t16 - 0.33e2 * t6 + 0.10e2;
	return u_xx + u_yy - f;
}

