#include <stdlib.h>		// malloc, calloc, realloc, free
#include <iostream>		// cout, endl
#include <math.h>		// sqrt, log, pow, exp
#include <ctime>

using namespace std;

int bookleaf_getgeom() {
	clock_t startTimer, endTimer;
	int nshape, nel, iel, ierr, repeatFactor, ii, irf, *ielnd;
	double ONEBYNINE, rms, *a1, *a2, *a3, *b1, *b2, *b3, *elvol, *cnwt, *elx, *ely;

	nshape = 4;
//	nel = 2250; //2250 = noh_small, 162000 = noh_med, 1296000 = noh_large, 2500 = sod_small, 250000 = sod_med, 25000000 = sod_large
	nel = 25000000; //2250 = noh_small, 162000 = noh_med, 1296000 = noh_large, 2500 = sod_small, 250000 = sod_med, 25000000 = sod_large
	ONEBYNINE=1.0/9.0;
	repeatFactor = 10;

	startTimer = clock();

	a1 = (double *) calloc (nel, sizeof(double));
	a2 = (double *) calloc (nel, sizeof(double));
	a3 = (double *) calloc (nel, sizeof(double));
	b1 = (double *) calloc (nel, sizeof(double));
	b2 = (double *) calloc (nel, sizeof(double));
	b3 = (double *) calloc (nel, sizeof(double));
	elvol = (double *) calloc (nel, sizeof(double));
	elx = (double *) calloc (nel*nshape, sizeof(double));
	ely = (double *) calloc (nel*nshape, sizeof(double));
	cnwt = (double *) calloc (nel*nshape, sizeof(double));
	ielnd = (int *) calloc (nel*nshape, sizeof(int));

	endTimer = clock();

	cout << "Allocate time taken: " << ((double)(endTimer - startTimer))/CLOCKS_PER_SEC << endl;

	startTimer = clock();


	for (int iel = 0; iel < nel; iel++) {
		for (int ii = 0; ii < nshape; ii++) {
			elx[iel*nshape+ii] = (double) (iel*ii*(7.0/3.0));
			ely[iel*nshape+ii] = (double) (iel*nshape+ii*(9.0/7.0));
		}
	}

	endTimer = clock();

	cout << "Initalize Values time taken: " << ((double)(endTimer - startTimer))/CLOCKS_PER_SEC << endl;

	startTimer = clock();


	for (int irf=0; irf < repeatFactor; irf++) {
//		! Calculate volume and iso-parametric terms
		for (int iel=0; iel < nel; iel++) {
			a1[iel]=0.25*(-elx[iel*nshape]+elx[iel*nshape+1]+elx[iel*nshape+2]-elx[iel*nshape+3]);
			a2[iel]=0.25*( elx[iel*nshape]-elx[iel*nshape+1]+elx[iel*nshape+2]-elx[iel*nshape+3]);
			a3[iel]=0.25*(-elx[iel*nshape]-elx[iel*nshape+1]+elx[iel*nshape+2]+elx[iel*nshape+3]);
			b1[iel]=0.25*(-ely[iel*nshape]+ely[iel*nshape+1]+ely[iel*nshape+2]-ely[iel*nshape+3]);
			b2[iel]=0.25*( ely[iel*nshape]-ely[iel*nshape+1]+ely[iel*nshape+2]-ely[iel*nshape+3]);
			b3[iel]=0.25*(-ely[iel*nshape]-ely[iel*nshape+1]+ely[iel*nshape+2]+ely[iel*nshape+3]);
			cnwt[iel*nshape]=ONEBYNINE*((3.0*b3[iel]-b2[iel])*(3.0*a1[iel]-a2[iel])-(3.0*a3[iel]-a2[iel])*(3.0*b1[iel]-b2[iel]));
			cnwt[iel*nshape+1]=ONEBYNINE*((3.0*b3[iel]+b2[iel])*(3.0*a1[iel]-a2[iel])-(3.0*a3[iel]+a2[iel])*(3.0*b1[iel]-b2[iel]));
			cnwt[iel*nshape+2]=ONEBYNINE*((3.0*b3[iel]+b2[iel])*(3.0*a1[iel]+a2[iel])-(3.0*a3[iel]+a2[iel])*(3.0*b1[iel]+b2[iel]));
			cnwt[iel*nshape+3]=ONEBYNINE*((3.0*b3[iel]-b2[iel])*(3.0*a1[iel]+a2[iel])-(3.0*a3[iel]-a2[iel])*(3.0*b1[iel]+b2[iel]));
			elvol[iel]=4.0*(a1[iel]*b3[iel]-a3[iel]*b1[iel]);
		}

		rms = 0.0;
		for (int ii=0; ii < nel; ii++) {
			rms = rms + elvol[ii];
		}

		std::cout << "repeatFactor:" << irf << "\trms*nel = " << rms << "\trms = " << sqrt(rms/(double)nel) << std::endl;
	}

	endTimer = clock();

	cout << "Computation time taken: " << ((double)(endTimer - startTimer))/CLOCKS_PER_SEC << endl;

	free(a1);
	free(a2);
	free(a3);
	free(b1);
	free(b2);
	free(b3);
	free(elvol);
	free(elx);
	free(ely);
	free(cnwt);
	free(ielnd);
}

/*
	The original data structure for this kernel was SOA:
		double * a1;
		double * a2; //...

	In order to put into SDLT, I first tried putting the arrays into a container, and using a container for all data:
		sdlt::soa1d_container<double> soaContainer(a1, nel);
		sdlt::soa1d_container<double> soaContainer(a2, nel); //ERROR: Redeclaration of soaContainer

	Therefore, to pass the data into SDLT, a structure has to be provided. As SDLT does not allow arrays or pointers in the structures,
	the elements have to be split up into 2 structures (those with "nel" elements, and those with "nel*nshape" elements). This is what
	has been done in this example.

	A 14% slowdown was seen using Intel 18.0 on a i5-3470 CPU. It was compiled using the follwoing line:
		icpc -O3 -std=c++11 bookleaf-getgeom-soa_sdlt.cpp 

	It should be noted that the value printed is consistent on runs. As dummy data is fed in (as we are conserned about that access
	patten and data layout affecting the wall clock time, rather than producing the correct result), the value produced is the same as
	the reference.
*/

#include <sdlt/sdlt.h>
#include <sdlt/soa1d_container.h>

struct bookleaf_data_1 {
	double a1,a2,a3,b1,b2,b3,elvol;
};

SDLT_PRIMITIVE(bookleaf_data_1, a1, a2, a3, b1, b2, b3, elvol);

struct bookleaf_data_2 {
	double elx,ely,cnwt,ielnd;
};

SDLT_PRIMITIVE(bookleaf_data_2, elx, ely, cnwt, ielnd);

int bookleaf_getgeom_sdlt() {
	//Variable declaration
	clock_t startTimer, endTimer;
	int nshape, nel, iel, ierr, repeatFactor, ii, irf, *ielnd;
	double ONEBYNINE, rms, *a1, *a2, *a3, *b1, *b2, *b3, *elvol, *cnwt, *elx, *ely;

	//Assign values
	nshape = 4;
//	nel = 2250; //2250 = noh_small, 162000 = noh_med, 1296000 = noh_large, 2500 = sod_small, 250000 = sod_med, 25000000 = sod_large
	nel = 25000000; //2250 = noh_small, 162000 = noh_med, 1296000 = noh_large, 2500 = sod_small, 250000 = sod_med, 25000000 = sod_large
	ONEBYNINE=1.0/9.0;
	repeatFactor = 10;

	//Allocate data
	startTimer = clock();

	sdlt::soa1d_container<bookleaf_data_1> bookleaf_container_1(nel);
	auto data_1 = bookleaf_container_1.access();

	sdlt::soa1d_container<bookleaf_data_2> bookleaf_container_2(nel*nshape);
	auto data_2 = bookleaf_container_2.access();

	endTimer = clock();

	std::cout << "Allocate time taken: " << ((double)(endTimer - startTimer))/CLOCKS_PER_SEC << std::endl;

	//Initialise Values
	startTimer = clock();

	for (int iel = 0; iel < nel; iel++) {
		for (int ii = 0; ii < nshape; ii++) {
			data_2[iel*nshape+ii].elx() = (double) (iel*ii*(7.0/3.0));
			data_2[iel*nshape+ii].ely() = (double) (iel*nshape+ii*(9.0/7.0));
		}
	}

	endTimer = clock();

	std::cout << "Initalize Values time taken: " << ((double)(endTimer - startTimer))/CLOCKS_PER_SEC << std::endl;

	//Compuation (access pattern correct)
	startTimer = clock();

	for (int irf=0; irf < repeatFactor; irf++) {
//		! Calculate volume and iso-parametric terms
		for (int iel=0; iel < nel; iel++) {
			data_1[iel].a1()=0.25*(-data_2[iel*nshape].elx()+data_2[iel*nshape+1].elx()+data_2[iel*nshape+2].elx()-data_2[iel*nshape+3].elx());
			data_1[iel].a2()=0.25*( data_2[iel*nshape].elx()-data_2[iel*nshape+1].elx()+data_2[iel*nshape+2].elx()-data_2[iel*nshape+3].elx());
			data_1[iel].a3()=0.25*(-data_2[iel*nshape].elx()-data_2[iel*nshape+1].elx()+data_2[iel*nshape+2].elx()+data_2[iel*nshape+3].elx());
			data_1[iel].b1()=0.25*(-data_2[iel*nshape].ely()+data_2[iel*nshape+1].ely()+data_2[iel*nshape+2].ely()-data_2[iel*nshape+3].ely());
			data_1[iel].b2()=0.25*( data_2[iel*nshape].ely()-data_2[iel*nshape+1].ely()+data_2[iel*nshape+2].ely()-data_2[iel*nshape+3].ely());
			data_1[iel].b3()=0.25*(-data_2[iel*nshape].ely()-data_2[iel*nshape+1].ely()+data_2[iel*nshape+2].ely()+data_2[iel*nshape+3].ely());
			data_2[iel*nshape].cnwt()=ONEBYNINE*((3.0*data_1[iel].b3()-data_1[iel].b2())*(3.0*data_1[iel].a1()-data_1[iel].a2())-(3.0*data_1[iel].a3()-data_1[iel].a2())*(3.0*data_1[iel].b1()-data_1[iel].b2()));
			data_2[iel*nshape+1].cnwt()=ONEBYNINE*((3.0*data_1[iel].b3()+data_1[iel].b2())*(3.0*data_1[iel].a1()-data_1[iel].a2())-(3.0*data_1[iel].a3()+data_1[iel].a2())*(3.0*data_1[iel].b1()-data_1[iel].b2()));
			data_2[iel*nshape+2].cnwt()=ONEBYNINE*((3.0*data_1[iel].b3()+data_1[iel].b2())*(3.0*data_1[iel].a1()+data_1[iel].a2())-(3.0*data_1[iel].a3()+data_1[iel].a2())*(3.0*data_1[iel].b1()+data_1[iel].b2()));
			data_2[iel*nshape+3].cnwt()=ONEBYNINE*((3.0*data_1[iel].b3()-data_1[iel].b2())*(3.0*data_1[iel].a1()+data_1[iel].a2())-(3.0*data_1[iel].a3()-data_1[iel].a2())*(3.0*data_1[iel].b1()+data_1[iel].b2()));
			data_1[iel].elvol()=4.0*(data_1[iel].a1()*data_1[iel].b3()-data_1[iel].a3()*data_1[iel].b1());
		}

		rms = 0.0;
		for (int ii=0; ii < nel; ii++) {
			rms = rms + data_1[ii].elvol();
		}

		std::cout << "repeatFactor:" << irf << "\trms*nel = " << rms << "\trms = " << sqrt(rms/(double)nel) << std::endl;
	}

	endTimer = clock();

	std::cout << "Computation time taken: " << ((double)(endTimer - startTimer))/CLOCKS_PER_SEC << std::endl;
}

int main () {
	clock_t startTimer;
	clock_t endTimer;
	startTimer = clock();
	bookleaf_getgeom();
	endTimer = clock();
	std::cout << "\\\\\\\\\\ Time taken: " << ((double)(endTimer - startTimer))/CLOCKS_PER_SEC << " /////" << std::endl << std::endl;
	startTimer = clock();
	bookleaf_getgeom_sdlt();
	endTimer = clock();
	std::cout << "\\\\\\\\\\ Time taken: " << ((double)(endTimer - startTimer))/CLOCKS_PER_SEC << " /////" << std::endl << std::endl;
	return 0;
}