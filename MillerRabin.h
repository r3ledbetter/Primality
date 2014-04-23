#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <vector>

#include <omp.h>

#include "mpreal.h"

using mpfr::mpreal;
using mpfr::mod;
using std::cout;
using std::endl;
using std::vector;

mpreal modularExponentiation(mpreal aVal, mpreal mVal, const mpreal& nVal);
mpreal getLargeInt(int digits);
mpreal getMaybePrime(int digits);

class MillerRabin {
	vector<mpreal> aValues;
	int iterations;
	int digits;
	int numLiars;

public:
	// for using small primes for a
	MillerRabin(const int numIter) {
		numLiars = 0;
		generateFirstNPrimes(numIter);
	}
	// for using random values for a
	MillerRabin(const int numIter, const int d, const double percentDigits, const mpreal& max) {
		digits = d;
		iterations = numIter;
		numLiars = 0;
		generateRandomA(percentDigits, max);
	}

	// maximum is first 100 primes
	void generateFirstNPrimes(const int numPrimes);
	// how many digits (as % of prime) in each a
	void generateRandomA(const double percentDigits, const mpreal& max);
	int getNumLiars() { return numLiars; }


	bool multThreadTest(const mpreal& possiblePrime);
	// single threaded millerRabin -> faster in most cases since we expect failure early for most numbers
	bool singleThreadTest(const mpreal& possiblePrime);
	// only run this test on KNOWN composites to try and find liars
	void liarTest(const mpreal& possiblePrime);
};

bool MillerRabin::singleThreadTest(const mpreal& possiblePrime) {
	// pull out all factors of 2 from possiblePrime-1
	if (possiblePrime == 1) return false;
	mpreal mOdd = possiblePrime - 1;
	int kFactorsOf2 = 0;
	while (mod(mOdd, 2) == 0) {
		++kFactorsOf2;
		mOdd = mOdd/2;
	}

	for (int i = 0; i < aValues.size(); ++i) {
		mpreal bVal = modularExponentiation(aValues[i], mOdd, possiblePrime);

		// initial check of computed bVal
		if ((bVal == 1) || (bVal == possiblePrime - 1)) {
			continue;
		}
		else if (kFactorsOf2 == 1) {
			return false;
		}
		// now check bVal after repetitive squaring until bIndex = kFactorsOf2
		for (int bIndex = 1; bIndex < kFactorsOf2; ++bIndex) {
			bVal = mod((bVal * bVal), possiblePrime);
			if (bVal == 1) {
				return false;
			}
			else if (bVal == possiblePrime - 1) {
				break;
			}
			else if (bIndex + 1 == kFactorsOf2) {
				return false;
			}
		}
	}
	// this means possiblePrime was even...
	if (kFactorsOf2 == 0) return false;
	
	return true;
}

void MillerRabin::liarTest(const mpreal& possiblePrime) {
	// pull out all factors of 2 from possiblePrime-1
	mpreal mOdd = possiblePrime - 1;
	int kFactorsOf2 = 0;
	while (mod(mOdd, 2) == 0) {
		++kFactorsOf2;
		mOdd = mOdd/2;
	}

	for (int i = 0; i < aValues.size(); ++i) {
		mpreal bVal = modularExponentiation(aValues[i], mOdd, possiblePrime);

		// initial check of computed bVal
		if ((bVal == 1) || (bVal == possiblePrime - 1)) {
			++numLiars;
			continue;
		}
		else if (kFactorsOf2 == 1) {
			continue; // this one isn't a liar
		}
		// now check bVal after repetitive squaring until bIndex = kFactorsOf2
		for (int bIndex = 1; bIndex < kFactorsOf2; ++bIndex) {
			bVal = mod((bVal * bVal), possiblePrime);
			if (bVal == 1) {
				break; // this one isn't a liar
			}
			else if (bVal == possiblePrime - 1) {
				++numLiars;
				break;
			}
			else if (bIndex + 1 == kFactorsOf2) {
				break; // we are done
			}
		}
	}
}

bool MillerRabin::multThreadTest(const mpreal& possiblePrime) {
	// pull out all factors of 2 from possiblePrime-1
	if (possiblePrime == 1) return false;
	mpreal mOdd = possiblePrime - 1;
	int kFactorsOf2 = 0;
	while (mod(mOdd, 2) == 0) {
		++kFactorsOf2;
		mOdd = mOdd/2;
	}

	volatile bool flag = false;
	bool returnVal = true;
	
	#pragma omp parallel for shared(flag)
	for (int i = 0; i < aValues.size(); ++i) {
		if (flag) continue;
		mpreal bVal = modularExponentiation(aValues[i], mOdd, possiblePrime);

		// initial check of computed bVal
		if ((bVal == 1) || (bVal == possiblePrime - 1)) {
			continue;
		}
		else if (kFactorsOf2 == 1) {
			returnVal = false;
			flag = true;
		}
		else {
			// now check bVal after repetitive squaring until bIndex = kFactorsOf2
			for (int bIndex = 1; bIndex < kFactorsOf2; ++bIndex) {
				bVal = mod((bVal * bVal), possiblePrime);
				if (bVal == 1) {
					returnVal = false;
					flag = true;
					break;
				}
				else if (bVal == possiblePrime - 1) {
					break;
				}
				else if (bIndex + 1 == kFactorsOf2) {
					returnVal = false;
					flag = true;
					break;
				}
			}
		}
	}
	// this means possiblePrime was even...
	if (kFactorsOf2 == 0) returnVal = false;
	
	return returnVal;
}

void MillerRabin::generateFirstNPrimes(const int numPrimes) {
	int bound = 0;
	if (numPrimes <= 10) bound = 30;
	else if (numPrimes <= 20) bound = 72;
	else if (numPrimes <= 30) bound = 114;
	else if (numPrimes <= 40) bound = 174;
	else bound = 542;

	vector<mpreal> smallPrimes;
	vector<int> eratos;
	eratos.push_back(0);
	eratos.push_back(0);
	eratos.push_back(2);
	for (int i = 3; i < bound; i += 2) {
		eratos.push_back(i);
		eratos.push_back(0);
	}
	for (int p = 0; p < bound; ++p) {
		if (eratos[p] != 0 && smallPrimes.size() < numPrimes) {
			p = eratos[p];
			smallPrimes.push_back(p);
			for (int delIndex = 2*p; delIndex < bound; delIndex +=p) {
				if (eratos[delIndex] != 0) {
					eratos[delIndex] = 0;
				}
			}
		}
	}
	aValues = smallPrimes;
}

void MillerRabin::generateRandomA(const double percentDigits, const mpreal& max) {
	int digitsInA = percentDigits * digits;
	mpreal decimalVal;
	mpreal largeInt;
	vector<mpreal> randomA;
	for (int i = 0; i < iterations; ++i) {
		decimalVal = mpfr::random(0);
		for (int j = 0; j < digitsInA; ++j) {
			decimalVal *= 10;
		}
		mpfr::modf(decimalVal, largeInt);
		// must be between [2, max-1]
		if (largeInt > 1 && largeInt < max) {
			randomA.push_back(largeInt);
		}
		else {
			--i;
		}
	}
	aValues = randomA;
}

// compute bVal = aVal^mVal (mod nVal)
mpreal modularExponentiation(mpreal aVal, mpreal mVal, const mpreal& nVal) {
	mpreal bVal = 1;
	while (mVal > 0) {
		if (mod(mVal, 2) == 0) {
			mVal = mVal/2;
			aVal = mod((aVal * aVal), nVal);
		}
		else {
			--mVal;
			bVal = mod((bVal * aVal), nVal);
		}
	}
	return bVal;
}

mpreal getLargeInt(int digits) {
	mpreal decimalVal = mpfr::random(0);
	for (int i = 0; i < digits; ++i) {
		decimalVal *= 10;
	}
	mpreal largeInt;
	mpfr::modf(decimalVal, largeInt);
	return largeInt;
}

// return a number with n digits that doesn't have small factors
mpreal getMaybePrime(int digits) {
	mpreal maybePrime;

	int smallFactors[10] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29};
	int done = false;
	while (!done) {
		maybePrime = getLargeInt(digits);
		done = true;
		for (int i = 0; i < 10; ++i) {
			if (mod(maybePrime, smallFactors[i]) == 0) {
				done = false;
				break;
			}
		}
	}
	return maybePrime;
}
