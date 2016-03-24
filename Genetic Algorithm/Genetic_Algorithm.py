import random
import dna_translator
import csv
from matplotlib import pyplot as plt
import numpy as np

"""
PREREQUISITE:
Create a directory 'GA Results Data' within main GA directory
So that your csv data from your GA runs can be stored.
"""

"""
~~~ GLOBAL VARIABLES ~~~
"""
# I chose bridge helix of CRISPR-associated endonuclease Cas9 because it
# is a small domain and CRISPR is a well-known bioligical system
OPTIMAL = 'EATRLKRTARRRYTRRKNRICYLQEIFSNEMAKVD'
OPTIMAL_AS_DNA = 'gaagcgacccgcctgaaacgcaccgcgcgccgccgctatacccgccgcaaaaaccgcatttgctatctgcaggaaatttttagcaacgaaatggcgaaagtggat'
DNA_SIZE = len(OPTIMAL) * 3
POP_SIZE = 6
GENERATIONS = 20000
MUTATION_CHANCE = 1.0/100.0
chromosome = 'a' + 't' + 'c' + 'g'


"""
~~~ GRAPH FUNCTION ~~~
"""
def plot_results(results_data):

	dat = np.array(results_data)
	plt.scatter(dat[:,0],dat[:,1],color="b")
	plt.grid()
	plt.title('Genetic Algorithm: Generation vs. Fitness')
	plt.ylabel('Fitness')
	plt.xlabel('Generation')
	plt.show()


"""
~~~ SUPPLEMENTARY FUNCTIONS ~~~
"""

# Creates a random string of nucleotides aka a single individual.
def random_individual():
    return ''.join((random.choice(chromosome)) for x in range(DNA_SIZE))

# Return a list of POP_SIZE individuals, each randomly generated via iterating
# DNA_SIZE times from the strings created from individual.
def random_population():
    pop_list = []
    for i in xrange(POP_SIZE):
        dna = random_individual()
        pop_list.append(dna)
    return pop_list


"""
~~~ GENETIC ALGORITHM FUNCTIONS ~~~
"""

# Compare the letters side by side. If they mismatch add 1 to the fitness value.
# The lower the number the better.
def fitness(dna):
    fitness = 0
    for i in xrange(DNA_SIZE):
        if dna[i] == OPTIMAL_AS_DNA[i]:
            fitness += 0
        elif dna[i] != OPTIMAL_AS_DNA[i]:
            fitness += 1
    return fitness

# Cross over a pair of DNA which would yield two offspring.
# Slices both dna1 and dna2 into two parts at a random index within their
# length and merges them. Both keep their initial sublist up to the crossover
# index, but their ends are swapped.
def crossover(dna1, dna2):
    position = random.randint(0, DNA_SIZE)
    return [dna1[:position]+dna2[position:], dna2[:position]+dna1[position:]]

# For every nucleotide in the DNA, there is a MUTATION_CHANCE chance that it will
# be swapped with a different nucleotide. This promotes diversity in the population.
def mutate(dna):
    dna_out = ""
    for i in xrange(DNA_SIZE):
        if random.random() > MUTATION_CHANCE:
            dna_out += dna[i]
        else:
            dna_out += random.choice(chromosome)
    return dna_out


"""
~~~ MAIN DRIVER OF GENETIC ALGORITHM ~~~
"""

if __name__ == "__main__":

	population = random_population()
	print "Original Population: ", population

	# This list for writing csv files.
	results_data = []
	for generation in xrange(GENERATIONS):

		# Make a sorted list of individuals of the population and their fitness in ascending order.
		sorted_pop_w_fitness = sorted([[indiv, fitness(indiv)] for indiv in population], key=lambda pair: pair[1])
		# print "List: sorted_pop_w_fitness =", sorted_pop_w_fitness

		sorted_individuals = []
		# Create list of sorted individuals from 'sorted_pop_w_fitness' without the fitness value.
		for x in xrange(POP_SIZE):
			sorted_individuals.append(sorted_pop_w_fitness[x][0])
		# print "List: sorted_individuals =", sorted_individuals

		"""
		Select, Crossover, Mutate
		"""
		# Selection
		# rtype: list(list(string,string))
		coupled_dna = []
		for a in xrange(0,POP_SIZE-1,2):
			coupled_dna.append([sorted_individuals[a],sorted_individuals[a+1]])
		# print "List: coupled_dna =", coupled_dna

		# Crossover
		# rtype: list(function(string,string))
		crossover_dna = []
		for a in xrange(POP_SIZE/2):
			crossover_dna.append(crossover(coupled_dna[a][0],coupled_dna[a][1]))
		# print "List: crossover_dna =", crossover_dna

		# Mutate
		flattened_list = [item for sublist in crossover_dna for item in sublist]
		for a in xrange(len(flattened_list)):
			flattened_list[a] = mutate(flattened_list[a])
			mutated_dna = flattened_list
		# print "List: mutated_dna =", mutated_dna

		"""
		Evaluate new fitness
		"""
		# Create a list of mutated individuals and their fitness
		post_sorted_pop_w_fitness = sorted([[indiv, fitness(indiv)] for indiv in mutated_dna], key=lambda pair: pair[1])
		# print "List: post_sorted_pop_w_fitness = ", post_sorted_pop_w_fitness

		post_sorted_individuals = []
		# Create list of sorted individuals from 'post_sorted_pop_w_fitness' without the fitness value.
		for b in xrange(POP_SIZE):
			post_sorted_individuals.append(post_sorted_pop_w_fitness[b][0])
		# print "List: post_sorted_individuals =", post_sorted_individuals

		for b in xrange(POP_SIZE):
			population[b] = post_sorted_individuals[0]
		# print population

		fittest_individual = population[0]
		fittest_protein = dna_translator.dna_to_protein(fittest_individual)

		results_data.append([generation,post_sorted_pop_w_fitness[0][1]])

		print "Generation %s ... Most Fit Individual: '%s'" % (str(generation).zfill(4), fittest_protein)
		print "                                Fitness: %i" % (post_sorted_pop_w_fitness[0][1])
		print ""

		if fitness(fittest_individual) == 0:
			break

	"""
	Display the final fittest individual
	"""
	print "FITTEST INDIVIDUAL: '%s'" % (fittest_protein)
	print "           FITNESS: %i" % (post_sorted_pop_w_fitness[0][1])
	print "        GENERATION: %i" % (generation)


	plot_results(results_data)

	"""
	Write a CSV file in directory 'GA Results Data' (located within main GA directory)
	"""
	fileptr = open('GA Results Data/popsize6.csv','w')
	writer = csv.writer(fileptr)
	for line in results_data:
		writer.writerow(line)
	fileptr.close()

	exit(0)
