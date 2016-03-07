/*
Example of how to compare exp(-BIGNUMBER) and exp(-OTHERBIGNUMBER) by using mpfr
*/
#include<iostream>
#include<sstream>
#include<mpfr.h>
#include<math.h>
#include<cmath>        // std::abs

// To compile and run: g++ exampleMpfr.cpp -o file -lmpfr -lgmp; ./file

int expMpFr(mpfr_t expOpDividedByTwoSquared, mpfr_t op,int prec) {
//The limit is between 744261110 and 744261120;

//int mpfr_cmp_d(mpfr_t op1, double op2) // > 0 if op1 > op2
mpfr_set_default_prec(prec);
mpfr_t two,opDividedByTwo,expOpDividedByTwo;
mpfr_inits2 (prec,two,opDividedByTwo,expOpDividedByTwo,(mpfr_ptr) 0); // and initialize them with precision 256

mpfr_set_str(two,"2.0",10,MPFR_RNDN);

mpfr_div(opDividedByTwo, op, two, MPFR_RNDN);

mpfr_exp(expOpDividedByTwo,opDividedByTwo,MPFR_RNDN);

mpfr_pow(expOpDividedByTwoSquared,expOpDividedByTwo,two,MPFR_RNDN);

mpfr_clears(two,opDividedByTwo,expOpDividedByTwo,(mpfr_ptr) 0);

return 1;
}

int main (int argc, char **argv) {

    long double truc=100844261120;//0; The limit is between 744261110 and 744261120;
    long double autreTruc=100844261121;//01;
    int prec=256;
    mpfr_rnd_t rnd=MPFR_RNDU;
    
    // == Convert doubles to str ==
    std::ostringstream strs;
    strs.precision(20);
    strs << -truc; //-truc/2;
    std::string strstr = strs.str();
    const char* str1 = strstr.c_str();
    std::ostringstream strs2;
    strs2.precision(20);
    strs2 << -autreTruc; //-autreTruc/2;
    std::string strstr2 = strs2.str();
    const char* str2 = strstr2.c_str();
    std::cout << "str1:" << str1 << std::endl;
    std::cout << "str2:" << str2 << std::endl;
    std::cout << std::endl;
    // == end of convert to doubles to str ==


    double nan=1.0/0.0;
    double aa=1.0e200;
    double bb=-1.0e200;
    bool hum=std::abs(aa*bb) >= nan;
    std::cout << std::abs(aa*bb) << " >= nan:" << hum << std::endl;
    
    mpfr_set_default_prec(prec);
    mpfr_t a,b,c,d,c2,d2; // Declare four numbers
    mpfr_inits2 (prec,a,b,c,d,c2,d2,(mpfr_ptr) 0); // and initialize them with precision 256

    mpfr_set_str(a,str1,10,rnd); // Variable where we store, value as a char*, decimal base, rounding option
    mpfr_set_str(b,str2,10,rnd); // Variable where we store, value as a char*, decimal base, rounding option
    // mpfr_set_str(a,"15", 10,MPFR_RNDN); // Variable where we store, value as a char*, decimal base, rounding option

    mpfr_printf ("a : %.128RNf \n", a);
    mpfr_printf ("b : %.128RNf \n", b);

    mpfr_exp(c,a,rnd);
    mpfr_exp(d,b,rnd);

    expMpFr(c2,a,prec);
    expMpFr(d2,b,prec);

    std::cout << "c:" << std::endl;
    mpfr_out_str(stdout,10,prec,c,rnd); // stdout, decimal base, precision, rounding option
    std::cout << std::endl;
    std::cout << "d:" << std::endl;
    mpfr_out_str(stdout,10,prec,d,rnd); // stdout, decimal base, precision, rounding option
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "c2:" << std::endl;
    mpfr_out_str(stdout,10,prec,c2,rnd); // stdout, decimal base, precision, rounding option
    std::cout << std::endl;
    std::cout << "d2:" << std::endl;
    mpfr_out_str(stdout,10,prec,d2,rnd); // stdout, decimal base, precision, rounding option
    std::cout << std::endl;
    std::cout << std::endl;

    mpfr_printf ("c : %.256RNe \n", c);
    mpfr_printf ("d : %.256RNe \n", d);
    mpfr_printf ("c2 : %.256RNe \n", c2);
    mpfr_printf ("d2 : %.256RNe \n", d2);
    
    std::cout.precision(prec);
    if (mpfr_cmp(c,d) == 1)
      std::cout << "exp(-" << truc << ") > exp(-" << autreTruc << ")" << std::endl;
    else
      std::cout << "exp(-" << truc << ") <= exp(-" << autreTruc << ")" << std::endl;
    
    mpfr_clears(a,b,c,d,c2,d2,(mpfr_ptr) 0); // Clear variables

    return 0;
}
