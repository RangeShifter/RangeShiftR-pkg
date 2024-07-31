#ifndef ALLELEH
#define ALLELEH

class Allele {
	const int id;
	const float value;
	const float dominance;
	inline static int counter = 0; // track nb alleles to set unique ID for each allele
public:
	Allele(float alleleValue, float alleleDominance) : value(alleleValue), dominance(alleleDominance), id(counter) { ++counter; }
	~Allele() {}
	float getAlleleValue() const { return value; };
	float getDominanceCoef() const { return dominance; };
	int getId() const { return id; }
};
#endif