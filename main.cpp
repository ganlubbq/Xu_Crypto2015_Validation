/* This is intended for validating the Appendix C in 
   Zhang B, Xu C, Meier W. Fast Correlation Attacks over Extension Fields, Large-Unit Linear Approximation and Cryptanalysis of SNOW 2.0[M]
     //Advances in Cryptology--CRYPTO 2015. Springer Berlin Heidelberg, 2015: 643-662.
   It seems they are wrong. Thus, the improvement is not as significant as declared.
   The NTL library is used for finite field arithmetics. The BOOST library is also used.

   Author: Yuan Yao
   Mail: yaoyuan1216@gmail.com
   Date: 2015-11-08
*/

#include <BerlekampMassey.h>
#include <NTL/GF2EX.h>
#include <NTL/GF2XFactoring.h>
#include <NTL/GF2EXFactoring.h>
#include <boost/lexical_cast.hpp>
#include <array>
#include <iostream>
#include <random>

int main() {
	using boost::lexical_cast;

	// Setup SNOW
	auto const snow = lexical_cast<NTL::GF2X>("[1 0 0 1 0 1 0 1 1]");
	BOOST_ASSERT(NTL::IterIrredTest(snow));
	NTL::GF2E::init(snow);

	NTL::GF2EX const snow_extended = [] {
		auto const beta = lexical_cast<NTL::GF2E>("[0 1]");
		NTL::GF2EX result;
		SetCoeff(result, 0, power(beta, 239));
		SetCoeff(result, 1, power(beta, 48));
		SetCoeff(result, 2, power(beta, 245));
		SetCoeff(result, 3, power(beta, 23));
		SetCoeff(result, 4);
		return result;
	}();
	BOOST_ASSERT(NTL::DetIrredTest(snow_extended));

	// Generate a random sequence
	size_t constexpr LFSR = 16;
	std::array<NTL::GF2EX, LFSR * 3> State;

	SetSeed(NTL::conv<NTL::ZZ>(std::random_device{}()));
	for (size_t i = 0; i != LFSR; ++i) {
		State[i] = NTL::random_GF2EX(deg(snow_extended));
	}
	auto const
		alpha = lexical_cast<NTL::GF2EX>("[[] [1]]"),
		alpha_inverse = InvMod(alpha, snow_extended);
	for (size_t i = LFSR; i != State.size(); ++i) {
		State[i] = (alpha_inverse * State[i - 5] + State[i - 14] + alpha * State[i - 16]) % snow_extended;
	}

	// Define test functions
	auto const print = [](auto const& f) {
		std::cout << "===================================================\n";
		std::cout << "Coefficients: \n";
		for (auto const& element : f) {
			std::cout << element << "\n";
		}
		std::cout << "Length: " << boost::size(f) << "\n";
		std::cout << "===================================================\n";
	};
	auto const Test = [&State, &print](auto const& M) {
		print(BerlekampMassey(State,
			[&](auto& a, auto const& b) {a = (a + b) % M;},
			[&](auto& a, auto const& b) {a = (a - b) % M;},
			[&](auto const& a, auto const& b) {return MulMod(a, b, M);},
			[&](auto const& a) {return InvMod(a, M);}
		));
	};

	// Reconstruct the LFSR on the SNOW field
	std::cout << "With SNOW\n";
	Test(snow_extended);

	// Reconstruct the LFSR on the AES field
	// Note that the binary representation of the sequence is not changed as expected.
	std::cout << "With AES\n";
	auto const aes = lexical_cast<NTL::GF2X>("[1 1 0 1 1 0 0 0 1]");
	BOOST_ASSERT(NTL::IterIrredTest(aes));
	NTL::GF2E::init(aes);
	NTL::GF2EX const aes_extended = lexical_cast<NTL::GF2EX>(
		"[[1 1 1 0 0 1 1] [0 0 1 0 0 1 1] [1 0 1 1 0 0 0 1] [0 0 1 1 1 0 1] [0 0 0 1 0 0 1] [1]]"
		);
	BOOST_ASSERT(NTL::DetIrredTest(aes_extended));
	Test(aes_extended);

	return 0;
}