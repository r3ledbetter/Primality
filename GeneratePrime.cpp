#include "MillerRabin.h"
#include <time.h>

mpreal getPrime(int digits, int iterations);

int DIGITS = 1000;

int main() {

	srand(time(NULL));
	mpfr::random(time(NULL));

	double wallTime;
	int maxDigits = 2 * DIGITS + 1;
	mpreal::set_default_prec(mpfr::digits2bits(maxDigits));
	cout.precision(maxDigits);

	// Time computations
	double lastIterLength = 0;
	double sum = 0;

	mpreal myPrime;
	cout << "Finding prime with " << DIGITS << " digits.." << endl;
	for (int i = 1; true; ++i) {
		wallTime = omp_get_wtime();
		myPrime = getPrime(DIGITS, 10);
		lastIterLength = (omp_get_wtime() - wallTime);
		sum += lastIterLength;
		cout << "Prime found: \n" << myPrime << "\n" << endl;
		cout << "Average time after " << i << " iterations: " << sum/i << "." << endl;
		cout << "Press enter to find another " << DIGITS << " digit prime." << endl;
		getchar();
	}
}

mpreal getPrime(int digits, int iterations) {
	int sizeofNumbers = 100;
	if (digits <= 10) {
		sizeofNumbers = 50;
	}
	else if (digits <= 50) {
		sizeofNumbers = 250;
	}
	else if (digits <= 200) {
		sizeofNumbers = 400;
	}
	else if (digits <= 300) {
		sizeofNumbers = 600;
	}
	else {
		sizeofNumbers = 2500;
	}

	//generate array of numbers in range [random, random+sizeofNumbers)
	mpreal* numbers = new mpreal[sizeofNumbers];
	numbers[0] = getLargeInt(digits);

	for (int i = 0; i < sizeofNumbers- 1; ++i) {
		numbers[i + 1] = numbers[i] + 1;
	}

	//seive of eratosthenes to eliminate numbers w/ small factors
	int smallFactors[10] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29};
	for (int f_index = 0; f_index < 10; ++f_index) {
		for (int n_index = 0; n_index < sizeofNumbers; ++n_index) {
			//check if already replaced at n_index
			if (numbers[n_index] == 0) continue;

			// check if small factor divides number
			if (mod(numbers[n_index], smallFactors[f_index]) == 0) {
				int factor = smallFactors[f_index];
				//set numbers[]=0 at intervals of factor 
				for (int interval = 0; (n_index + (interval * factor)) < sizeofNumbers; ++interval) {
					numbers[n_index + (interval * factor)] = 0;
				}
				break;
			}
		}
	}

	MillerRabin mr(iterations);
	// determine if mult threads would improve performance
	if (digits >= 8) {
		volatile bool flag = false;
		mpreal returnVal = 0;
		// use mult threads to run MR test on remaining numbers
		#pragma omp parallel for shared(flag)
		for (int i = 0; i < sizeofNumbers; ++i) {
			if (flag) continue;
			if (numbers[i] != 0) {
				if (mr.singleThreadTest(numbers[i])) {
					flag = true;
					returnVal = numbers[i];
				}
			}
		}
		if (returnVal == 0) {
			//cout << "Trying again..." << endl;
			return getPrime(digits, iterations);
		}
		else return returnVal;
	}
	else {
		// use single thread to run MR test on remaining numbers
		for (int i = 0; i < sizeofNumbers; ++i) {
			if (numbers[i] != 0) {
				if (mr.singleThreadTest(numbers[i])) {
					return numbers[i];
				}
			}
		}
		//cout << "Trying again..." << endl;
		return getPrime(digits, iterations);
	}
}
