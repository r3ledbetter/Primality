#include "../GeneratePrime/GeneratePrime.cpp"
#include <time.h>
#include <fstream>

void liarsWithRandomsVsSmallPrimes();
vector<mpreal> generateComposites(int range, MillerRabin mrShort);
void randomVsSmallPrimes();
void liarsOnCarmichaels();
bool isCarmichael(const mpreal& num);

int DIGITS = 200;

using std::cin;
using std::ifstream;

int main() {
	srand(time(NULL));
	mpfr::random(time(NULL));

	int maxDigits = 2*DIGITS + 1;
	mpreal::set_default_prec(mpfr::digits2bits(maxDigits));
	cout.precision(DIGITS);

	//liarsWithRandomsVsSmallPrimes();
	//randomVsSmallPrimes();
	liarsOnCarmichaels();
}

void liarsOnCarmichaels() {
	mpreal tempCarm;
	vector<mpreal> carmichaels;
	ifstream carmFile("CarmichaelNumbers.txt");
	for (int i = 0; i < 885; ++i) {
		carmFile >> tempCarm;
		carmichaels.push_back(tempCarm);
		// Initially check and make sure all numbers are carmichael numbers
		//if (!isCarmichael(carmichaels[i])) cout << i << endl;
	}
	carmFile.close();

	vector<int> primeStats(116);
	vector<int> randomStats(116);
	vector<int> semiRandomStats(116);
	vector<int> numCarmAtDigits(116);

	for (int i = 0; i < carmichaels.size(); ++i) {
		// Get number of digits in carmichaels[i]
		int digits = 0;
		mpreal temp = carmichaels[i];
		while (temp > 1) {
			temp /= 10;
			++digits;
		}

		// Find liars
		MillerRabin mrPrimes(100);
		mrPrimes.liarTest(carmichaels[i]);

		MillerRabin mrRandoms(100, digits, 1.0, carmichaels[i]);
		mrRandoms.liarTest(carmichaels[i]);

		MillerRabin mrSemiRandoms(100, digits, 0.5, carmichaels[i]);
		mrSemiRandoms.liarTest(carmichaels[i]);

		primeStats[digits] += mrPrimes.getNumLiars();
		randomStats[digits] += mrRandoms.getNumLiars();
		semiRandomStats[digits] += mrSemiRandoms.getNumLiars();
		++numCarmAtDigits[digits];
	}
	
	for (int i = 0; i < primeStats.size(); ++i) {
		if (numCarmAtDigits[i] == 0) continue;
		cout << "Stats for " << i << " digits. Carmichael numbers used: " << numCarmAtDigits[i] << endl
			<< " Prime liars total: " << primeStats[i] << endl
			<< " Random liars total: " << randomStats[i] << endl
			<< " \"Random\" liars total: " << semiRandomStats[i] << endl;
	}
}

// Compare number of liars found using small primes, random numbers, and "random" half-digit numbers for a
void liarsWithRandomsVsSmallPrimes() {
	int count = 0;
	int primeLiars = 0;
	int randomLiars = 0;
	int semiRandomLiars = 0;

	// This gets passed to generateComposites() and is used to find composites
	// It needs to be declared outside parallel section to work properly
	// Only need one iteration to guarantee a composite
	MillerRabin mrShort(1);

	#pragma omp parallel shared(count, primeLiars, randomLiars, semiRandomLiars, mrShort)
	{
		while (count < 100000000) {
			vector<mpreal> composites = generateComposites(1000, mrShort);

		/******************************************/

			// testing with small primes for a
			MillerRabin mrPrimes(100);
			for (int i = 0; i < composites.size(); ++i) {
				mrPrimes.liarTest(composites[i]);
			}
			primeLiars += mrPrimes.getNumLiars();

		/******************************************/

			// testing with random numbers for a
			MillerRabin mrRandoms(100, DIGITS, 1.0, composites[composites.size() - 1]);
			for (int i = 0; i < composites.size(); ++i) {
				mrRandoms.liarTest(composites[i]);
			}
			randomLiars += mrRandoms.getNumLiars();
		/******************************************/

			// testing with "random" (half DIGITS) numbers for a
			MillerRabin mrSemiRandoms(100, DIGITS, 0.5, composites[composites.size() - 1]);
			for (int i = 0; i < composites.size(); ++i) {
				mrSemiRandoms.liarTest(composites[i]);
			}
			semiRandomLiars += mrSemiRandoms.getNumLiars();

		/******************************************/

			count += 100*composites.size();	// b/c we do 100 MR tests on each composite
			cout << "PL: " << primeLiars << ", RL: " << randomLiars << ", SRM: " << semiRandomLiars << " -> Total count: " << count << endl;
		}
	}
	cout.precision(10);
	cout << "(Lower percentage is better) " << DIGITS << endl;
	cout << primeLiars << ", Percent primeLiars: " << (double)primeLiars/count*100 << endl;
	cout << randomLiars << ", Percent randomLiars: " << (double)randomLiars/count*100 << endl;
	cout << semiRandomLiars << ", Percent semiRandomLiars: " << (double)semiRandomLiars/count*100 << endl;
}

// get numbers without small factors that are known composites
vector<mpreal> generateComposites(int range, MillerRabin mrShort) {
	vector<mpreal> composites;
	// need to make sure 1 isn't in the composites vector
	mpreal firstNumber = 1;
	while (firstNumber == 1) {
		firstNumber = getLargeInt(DIGITS);
	}
	composites.push_back(firstNumber);

	for (int i = 0; i < range - 1; ++i) {
		composites.push_back(composites[i] + 1);
	}
	//seive of eratosthenes to eliminate numbers w/ small factors
	int smallFactors[10] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29};
	for (int f_index = 0; f_index < 10; ++f_index) {
		for (int n_index = 0; n_index < range; ++n_index) {
			if (composites[n_index] == 0) continue;

			// check if small factor divides number
			if (mod(composites[n_index], smallFactors[f_index]) == 0) {
				int factor = smallFactors[f_index];
				//set composites[]=0 at intervals of factor 
				for (int interval = 0; (n_index + (interval * factor)) < range; ++interval) {
					composites[n_index + (interval * factor)] = 0;
				}
				break;
			}
		}
	}

	// now we only want numbers that are NOT prime
	vector<mpreal> returnVec;
	for (int i = 0; i < range; ++i) {
		if (composites[i] != 0) {
			if (!mrShort.singleThreadTest(composites[i])) {
				returnVec.push_back(composites[i]);
			}
		}
	}
	return returnVec;
}

// compare execution times of MR test using random a vs small prime a
// single threaded only
void randomVsSmallPrimes() {
	while (true) {
		int iterations;

		cout << "\nNumber of iterations: ";
		cin >> iterations;
		cout << "Number of digits: ";
		cin >> DIGITS;

		cout << "Finding a(n) " << DIGITS << " digit prime number... " << endl;
		mpreal prime = getPrime(DIGITS, 40);
		cout.precision(DIGITS);
		cout << prime << endl;
		cout.precision(6);

		MillerRabin mrPrimes(10);
		MillerRabin mrRandoms(10, DIGITS, 1.0, prime);
		MillerRabin mrSemiRandoms(10, DIGITS, 0.5, prime);

		double totalTime;
	/******************************************/
		cout << "Starting test with small primes..." << endl;
		clock_t start, end;
		start = clock();

		for (int i = 0; i < iterations; ++i) {
			mrPrimes.singleThreadTest(prime);
		}

		end = clock();
		totalTime = (double)(end-start)/CLOCKS_PER_SEC;
		cout << "Time with small primes: " << totalTime/iterations << endl;

	/******************************************/
		cout << "Starting test with random numbers..." << endl;
		start = clock();

		for (int i = 0; i < iterations; ++i) {
			mrSemiRandoms.singleThreadTest(prime);
		}

		end = clock();
		totalTime = (double)(end-start)/CLOCKS_PER_SEC;
		cout << "Time with \"random\" numbers: " << totalTime/iterations << endl;

	/******************************************/
		cout << "Starting test with semirandom numbers..." << endl;
		start = clock();

		for (int i = 0; i < iterations; ++i) {
			mrRandoms.singleThreadTest(prime);
		}

		end = clock();
		totalTime = (double)(end-start)/CLOCKS_PER_SEC;
		cout << "Time with random numbers: " << totalTime/iterations << endl;
	}
}

// Make sure we are using carmichael numbers
// i.e. make sure they pass Fermat test and fail Miller-Rabin
bool isCarmichael(const mpreal& num) {
	MillerRabin mr(100); // strong test
	if (mr.singleThreadTest(num)) {
		// Uh oh, it's prime
		return false;
	}
	if (modularExponentiation(2, num-1, num) != 1) {
		// Uh oh, it's not a charmichael number
		return false;
	}
	return true;
}
