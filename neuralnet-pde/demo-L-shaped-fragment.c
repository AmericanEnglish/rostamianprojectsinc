/* The cut-off function for demo-L-shaped.c is
 *     phi(x,y) = (1-x^2)*(1-y^2)*(x + y + sqrt(x^2+y^2)).
 * The derivatives of phi where calculated symbolically
 * in Maple.  The result was translated to C code also in Maple.
*/
static void my_phi(struct Neural_Net_PDE *nn, double x, double y)
{
	double t1 = x * x;
	double t2 = t1 - 0.1e1;
	double t4 = y * y;
	double t5 = t4 - 0.1e1;
	double t6 = t1 + t4;
	double t7 = sqrt(t6);
	double t8 = t7 + y + x;
	double t9 = t8 * t5;
	double t11 = t5 * x;
	double t14 = t5 * t2;
	double t15 = 0.1e1 / t7;
	double t17 = t15 * x + 0.1e1;
	double t21 = y * t2;
	double t25 = t15 * y + 0.1e1;
	double t33 = 0.1e1 / t7 / t6;

	nn->phi[0][0] = t9 * t2;
	nn->phi[1][0] = 0.2e1 * t8 * t11 + t17 * t14;
	nn->phi[0][1] = t25 * t14 + 0.2e1 * t8 * t21;
	nn->phi[2][0] = 0.2e1 * t9 + 0.4e1 * t17 * t11 + (-t33 * t1 + t15) * t14;
	nn->phi[1][1] = -y * t33 * x * t14 + 0.4e1 * t8 * y * x + 0.2e1 * t25 * t11 + 0.2e1 * t17 * t21;
	nn->phi[0][2] = 0.2e1 * t8 * t2 + 0.4e1 * t25 * t21 + (-t33 * t4 + t15) * t14;
}

/* For the exact solution we pick u(x,t) = x*phi(x,y),
 * and then we "reverse engineer" the PDE by picking
 * f = u_xx + u_yy.  The calculation of f was done
 * symbolically in Maple.
*/
static double my_pde(double x, double y, double u,
		double u_x, double u_y,
		double u_xx, double u_xy, double u_yy)
{
	double t1 = x * x;
	double t2 = t1 * t1;
	double t6 = y * y;
	double t17 = sqrt(t1 + t6);
	double t22 = t6 * t6;
	double f = 0.2e1 / t17 * (t17 * (t2 + 0.3e1 * y * t1 * x + t1 * (0.6e1 * t6 - 0.7e1) + x * (0.3e1 * t6 * y - 0.6e1 * y) - t6 + 0.1e1) + x * (t2 + t1 * (0.19e2 / 0.2e1 * t6 - 0.15e2 / 0.2e1) + 0.3e1 * t22 - 0.15e2 / 0.2e1 * t6 + 0.3e1 / 0.2e1));
	return u_xx + u_yy - f;
}

