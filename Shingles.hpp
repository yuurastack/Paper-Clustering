#ifndef SHINGLES_HPP_INCLUDED
#define SHINGLES_HPP_INCLUDED

#include <vector>
#include <cstddef>      // NULL, std::size_t
#include <string>
#include <cstdlib>      // std::rand, std::srand
#include <ctime>        // std::time

using namespace std;

class Shingles {
public:
    typedef unsigned int Signature;

    Shingles();
    Shingles(int numhashes);

    template<typename ContainerT>
    Signature compute(const ContainerT&) const;

private:
    static const unsigned long bigPrimeNumber = 0x7FFFFFFF; // Same as 2^31 - 1
    unsigned int A;     // A, as in a X + b mod bigPrimeNumber
    unsigned int B;     // B, as in a X + b mod bigPrimeNumber

    vector<unsigned int> A_;  // A_, as in a X + b mod prime_
    vector<unsigned int> B_;  // B_, as in a X + b mod prime_
    vector<unsigned int> MIN_;// current MIN_HASH

    int numberSignatures;


};



inline
Shingles::Shingles() {
    std::srand(std::time(NULL));
    A = (std::rand() % bigPrimeNumber) + 1;
    B = (std::rand() % bigPrimeNumber) + 1;
    numberSignatures = 1;
}

Shingles::Shingles(int nhashes){
    numberSignatures = nhashes;
    srand(time(NULL));

    for (int i = 0; i < numberSignatures; i++){

      A_.push_back( rand() % bigPrimeNumber  + 1);
      B_.push_back( rand() % bigPrimeNumber + 1);
      MIN_.push_back( bigPrimeNumber + 1);
    }

}

template<typename ContainerT>
Shingles::Signature
Shingles::compute(const ContainerT& seqs) const {

    std::hash<typename ContainerT::value_type> hashing;
    vector<Signature> multiSignatures;
    Signature minShingleHash = bigPrimeNumber;

    for(int i=0; i<numberSignatures; i++){
        minShingleHash = bigPrimeNumber;
        for (const auto& val : seqs) {
            std::size_t hashedWord = hashing(val);
            std::size_t shingleID = hashedWord;
            //Signature shingleHash = (((unsigned long) A * (unsigned long) shingleID) + B) % bigPrimeNumber;
            Signature shingleHash = (((unsigned long) A_[i] * (unsigned long) shingleID) + B_[i]) % bigPrimeNumber;

            if (shingleHash < minShingleHash) {
                minShingleHash = shingleHash;
            }
        }
	multiSignatures.push_back(minShingleHash);	
    }
    for(unsigned int i=0; i<multiSignatures.size(); i++){
	//cout<<multiSignatures[i]<<" ";
    }	

    return minShingleHash;
}


#endif  // SHINGLES_HPP_INCLUDED
