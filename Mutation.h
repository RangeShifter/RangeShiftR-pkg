#ifndef MUTATIONH
#define MUTATIONH

class Mutation {
	const int id;
	const float selCoef;
	const float dominance;
	inline static int counter = 0;
public:
	Mutation(float selCoefA, float dominanceA) : selCoef(selCoefA), dominance(dominanceA), id(counter) { ++counter; }
	~Mutation() {}
	float getSelectionCoef() const { return selCoef; };
	float getDominanceCoef() const { return dominance; };
	float getId() const { return id; }
};
#endif