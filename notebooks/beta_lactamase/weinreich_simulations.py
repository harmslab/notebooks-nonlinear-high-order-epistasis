# This is a major hack...
#
# ------------------------------------------------
# Global variables
# ------------------------------------------------

MUTATIONS = {
            0:"g4205a",
            1:"A42G",
            2:"E104K",
            3:"M182T",
            4:"G238S",
            }

PATH_TO_DATA = "../../datasets/"

# ------------------------------------------------
# Imports
# ------------------------------------------------

import matplotlib.pyplot as plt
from matplotlib.cm import gray, spring 
import numpy as np
from collections import Counter
import networkx as nx

# ------------------------------------------------
# Initial loading of data.
# ------------------------------------------------

def load_weinrich_data(filename, raw=False):
    """ loads data from `filename` and returns genotypes-phenotypes """
    data = np.loadtxt(PATH_TO_DATA + filename, dtype="S20", delimiter=' ', ndmin=2)
    sequences = data[:,0].astype(str)
    fitness = np.array(data[:,1].astype(str), dtype=float)
    rank = np.array(data[:,4].astype(str), dtype=int)

    errors = np.ones(len(sequences), dtype=float) *  1.0/np.sqrt(3)
    for i in range(len(data)):
        vals = np.array(data[i, 1:4].astype(str), dtype=float)
        vals2 = list()
        for j in range(len(vals)):
            vals2.append(float(vals[j]))
        if raw:
            fitness[i] = sum(vals2)/3.0
        else:
            fitness[i] = math.log(sum(vals2)/3.0, math.sqrt(2)) # Log base sqrt(2)
  
    return sequences, fitness, errors, rank

sequences, fitnesses, errors, ranks = load_weinrich_data("weinrich_lactamase.txt")
title = "Weinrich Lactamase data"

# ------------------------------------------------
# Weinrich's probability distribution
# ------------------------------------------------

def probability_func(fitness1, fitness2):
    """ A probability function for Weinreich's data. """
    sij = (fitness2 - fitness1)/fitness1
    if sij < 0:
        sij = 0
    fixation = 1 - np.exp(-2*sij)
    return fixation

# ------------------------------------------------
# Methods for drawing trajectories on networks
# ------------------------------------------------

def edge_arrows(pos, edges):
    """ Maker a list of edge arrows. """
    arrows = list()
    for e in edges:
        arrows.append((pos[e[0]][0], pos[e[0]][1], 
                       .8*(pos[e[1]][0]-pos[e[0]][0]), 
                       .8*(pos[e[1]][1]- pos[e[0]][1]), 
                       edges[e]))
    return arrows

def edge_weight(traj):
    """ Count the number of times each edge is visited. """
    edge_counter = dict()
    for t in traj:
        sequences = t.split(",")
        edges = [(sequences[i-1],sequences[i]) for i in range(1,len(sequences))]
        for e in edges:
            if e in edge_counter:
                edge_counter[e] += traj[t]
            else:
                edge_counter[e] = traj[t]
    return edge_counter

def draw_space(G, pos=None):
    """ Draw the trajectories on Graph. """
    fig = plt.figure(figsize=[7,7])
    if pos is None:
        pos = nx.spring_layout(G, iterations=150)
     
    # Draw network
    colors = list()
    for n in G.nodes():
        colors.append(G.node[n]["phenotype"])
                      
    nx.draw(G,pos, alpha=.8, 
            cmap=gray, 
            node_color=colors, 
            node_size=400, 
            arrows=False, 
            with_labels=True, 
            width=0.5,
            vmin = 0.94,
            vmax = 1.2,
           )
    return pos, fig

def draw_traj(G, traj, pos=None):
    """ Draw the trajectories on Graph. """
    fig = plt.figure(figsize=[7,7])
    if pos is None:
        pos = nx.spring_layout(G, iterations=150)
     
    # Draw network
    colors = list()
    for n in G.nodes():
        colors.append(G.node[n]["phenotype"])
                      
    nx.draw(G,pos, alpha=.8, 
            cmap=gray, 
            node_color=colors, 
            node_size=400, 
            arrows=False, 
            with_labels=True, 
            width=2,
            vmin = 0.94,
            vmax = 1.2,
           )
    
    # Draw arrows
    edges = edge_weight(traj)
    arrows = edge_arrows(pos, edges)
    for a in arrows:
        plt.arrow(a[0], a[1], a[2], a[3], alpha=0.6, width=0.005*np.log(a[4]), head_width=0.05, head_length=0.05, fc='b', ec='k')
    return pos, fig

# ------------------------------------------------
# Fitness Class for Weinreich's simulation
# ------------------------------------------------

class WeinreichFitnesses(object):
    
    def __init__(self, sequences, ranks):
        
        self._sequences = sequences
        self._ranks = ranks
        
    @property
    def sequences(self):
        return self._sequences

    @property
    def ranks(self):
        return self._ranks

    @property
    def samples(self):
        """ Return fitness samples. """
        return self._samples

    @property
    def fitnesses(self):
        """ Return the mean fitnesses of each sequence."""
        return self._fitnesses

    @property
    def errors(self):
        """ Return the standard error for each sequence. """
        return self._errors
    
    @property
    def G(self):
        """ Return the networkx graph of this sequence space"""
        return self._G

    def generate_fitnesses(self, n_samples):
        """ Generate fitnesses from distribution. """
        self._samples = self._sample_fitness(0.1, n_samples)
        fits, errs = self._fitness_stats(self.samples)
        fitnesses = np.empty(len(self.sequences), dtype=float)
        errors = np.empty(len(self.sequences), dtype=float)
        for i in range(len(self.sequences)):
            # Find fitness rank and grab the value of that fitness element index
            fitnesses[i] = fits[self.ranks[i]-1] 
            errors[i] = errs[self.ranks[i]-1]
        self._fitnesses = fitnesses
        self._errors = errors

    # -------------------------------------
    # Hidden Methods
    # -------------------------------------

    def _fitness_dist(self, scale):
        """ Return 14 fitness values with average value equal to scale. """
        return np.sort(1 + np.random.exponential(scale=scale, size=14))[::-1]

    def _sample_fitness(self, scale, n_samples):
        """ Sample the fitness distribution n times."""
        samples = np.empty((14,n_samples))
        for i in range(n_samples):
            samples[:,i] = self._fitness_dist(scale)
        return samples

    def _fitness_stats(self, samples):
        """ Get mean fitnesses from a sample and their standard deviation. """
        return np.mean(samples, axis=1), np.std(samples, axis=1)/np.sqrt(samples.shape[1])

# -----------------------------------------------------
# Generic Class for Brute force Monte Carlo simulations
# of evolutionary trajectories in sequence space.
# -----------------------------------------------------

class SequenceSpace(object):
    
    def __init__(self, sequences, fitnesses, prob_func, mutations={}):
        """ Build a sequence space object on which trajectories can be calculated
            
            Args:
            ----
            sequences: array
                list of sequences
            fitnesses: array
                list of fitnesses for sequences
            prob_func: callable
                function for calculating transition probabilities from sequence
            mutations: dict
                mapping mutation index to string representation
        """
        self._sequences = sequences
        self._fitnesses = fitnesses
        self._indices = np.arange(len(sequences))
        self._mutations = mutations
        self.prob_func = prob_func
        self._build_graph()

        
    @property
    def sequences(self):
        return self._sequences
    
    @property
    def fitnesses(self):
        return self._fitnesses
    
    @property
    def indices(self):
        return self._indices
    
    @property
    def mutations(self):
        return self._mutations
    
    @property
    def G(self):
        """ Return the networkx graph of this sequence space"""
        return self._G    
    
    @property
    def i2seq(self):
        return dict(zip(self.indices, self.sequences))
    
    @property
    def seq2fit(self):
        return dict(zip(self.sequences, self.fitnesses))
   
    # -------------------------------------
    # Hidden Methods
    # -------------------------------------
    
    
    def _binary_neighbors(self, genotype):
        """ Returns binary genotypes that differ by a site in binary sequence space."""
        dim = len(genotype)
        chars = list(genotype)

        neighbors = list()
        for c in range(0,dim):
            nb = list(genotype)
            # Create a neighbor
            if chars[c] == '0':
                nb[c] = '1'
            else:
                nb[c] = '0'
            seq = "".join(nb)

            neighbors.append(seq)
        return neighbors

    def _compare_seq(self, seq1, seq2):
        """ return index of difference in sequence. """
        for i in range(len(seq1)):
            if seq1[i] != seq2[i]:
                index = i
        return index
        
    def _draw_edges(self, genotypes):
        """ Draw edges between binary genotypes. """
        edges = list()
        for g in genotypes:
            edges += [(g, neighbor) for neighbor in self._binary_neighbors(g)]
        return edges
    
    def _get_matrix(self):
        """ Return the transition matrix for G. """
        return nx.attr_matrix(self.G, edge_attr='probability', normalized=True, rc_order=self.sequences)
        
    def _probability_to_distribution(self, probs):
        """ Return the a list of probabilities as partitioned sections from 0 to 1. """
        return np.array([sum(probs[:p+1]) for p in range(len(probs))])
    
    def _to_sequence(self, trajectory):
        """ Convert a trajectory set of indices to binary sequences."""
        return [self.i2seq[t] for t in trajectory]
    
    def _to_mutations(self, trajectory):
        """ Convert a trajectory set of indices to mutations."""
        return [self.mutations[self._compare_seq(trajectory[i-1],trajectory[i])] for i in range(1,len(trajectory))]

    def trajectory(self):
        """ Traverse G from state=0 to state=len(G) based on transition matrix. Return trajectory."""
        matrix = self._get_matrix()
        max_iter = 1000
        counter = 0
        position = 0
        end = len(matrix)-1
        finished = False
        trajectory = (0,)
        while finished is False and counter < max_iter:
            # Separate transition probabilities into separate probabilty chucks to sample
            options = np.array(matrix[trajectory[-1]])[0]
            regions = self._probability_to_distribution(options)
            # Random number between 0 and 1
            rando = np.random.rand()
            new_pos = -1
            check = -1
            while new_pos == -1:
                check += 1
                if rando < regions[check]:
                    new_pos = check
                    trajectory += (new_pos,)
            if new_pos == end:
                finished = True
            counter +=1

        if counter == max_iter:
            raise Exception("Hit max iterations.")

        return trajectory
    
    def _build_graph(self):
        self._G = nx.DiGraph()
        edges = self._draw_edges(self.sequences)
        full_edges = list()
        gpm = self.seq2fit
        
        # Build list of network edges with their Prob(i-->j)
        for e in edges:
            seq1, seq2 = e[0], e[1]
            fitness1 = gpm[seq1]
            fitness2 = gpm[seq2]
            probability = self.prob_func(fitness1, fitness2)
            full_edges.append((seq1, seq2, {"probability":probability}))
        self._G.add_edges_from(full_edges)
        nx.set_node_attributes(self._G, 'phenotype', gpm)
    
    def enumerate_trajectories(self, n, binary=True):
        """ Enumerates n number of Monte Carlo trajectories and counts trajectories. """
        trajectories = list()
        for i in range(n):
            traj = self.trajectory()
            elements = self._to_sequence(traj)
            if binary is False:
                elements = self._to_mutations(elements)
            string = ",".join(elements)
            trajectories.append(string)
        return Counter(trajectories)