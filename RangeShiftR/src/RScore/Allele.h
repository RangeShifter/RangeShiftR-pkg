#ifndef ALLELEH
#define ALLELEH

class Allele {
	const float value;
	const float dominance;
public:
	Allele(float alleleValue, float alleleDominance) : value(alleleValue), dominance(alleleDominance) { }
	~Allele() {}
	float getAlleleValue() const { return value; };
	float getDominanceCoef() const { return dominance; };
};
#endif