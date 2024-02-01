#ifndef ALLELEH
#define ALLELEH

class Allele {
	const int id;
	const float value;
	const float dominance;
	inline static int counter = 0;
public:
	Allele(float alleleValue, float alleleDominance) : value(alleleValue), dominance(alleleDominance), id(counter) { ++counter; }
	~Allele() {}
	float getAlleleValue() const { return value; };
	float getDominanceCoef() const { return dominance; };
	float getId() const { return id; }
};
#endif