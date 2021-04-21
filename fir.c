double* fir;
double* sig;

int len_fir;


double convolute(double* x, double* h,int l) {
	double y;
	for (int i=l-1;i>=0;i-- ){
		y=y+(x[i]*h[i]);
	}
	y=y/l;
	return y
}
