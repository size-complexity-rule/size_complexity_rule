from random import randrange, choice, random

import sys, traceback

from math import exp, pow, log

from sys import argv



k          = int(argv[1])
k_rep      = float(argv[2])
k_dea      = float(argv[3])
stickiness = float(argv[4])
C          = float(argv[5])
p          = int(argv[8])

initial_population = int(argv[6])
t_max              = int(argv[7])


class Cell(object):
    def __init__(self, k, parent):
        self.k = k
        self.kids = []
        self.parent = parent
        self.total = 1

    def reproduce_random_member(self):
        found = False
        while not found:
            selected = randrange(self.total_cells())
            found = self.reproduce(selected)

    def reproduce(self, index=0):
        cell = self.search(index)
        if len(cell) < cell.k:
            cell.kids.append(Cell(cell.k, cell))
            cell.update_total_cells()
            return True
        else:
            return False

    def reorganize(self, l):
        if type(self.parent) is Population:
            self.parent.transfer(self, l)
        else:
            self.parent.reorganize(l)

    def root(self):
        if type(self.parent) is Population:
            return self
        else:
            return self.parent.root()

    def search(self, index):
        if index < 0:
            return "Failed search: negative index"

        elif index == 0:
            return self

        else:
            index -= 1
            for kid in self.kids:
                total_cells = kid.total_cells()
                if index < total_cells:
                    return kid.search(index)
                else:
                    index -= total_cells
            return None



    def kill_random_member(self):
        selected = randrange(self.total_cells())
        return self.kill(selected)

    def kill(self, index):
        cell = self.search(index)

        lst = []
        for kid in cell.kids:
            kid.parent = None
            lst.append(kid)

        if type(cell.parent) is not Population:
            cell.parent.kids.remove(cell)
            cell.update_total_cells()
            lst.append(cell.root())

        return lst

    def total_cells(self):
        return self.total

    def update_total_cells(self):
        self.total = 1 + sum([kid.total_cells() for kid in self.kids])
        if type(self.parent) is Cell:
            self.parent.update_total_cells()

    def __len__(self):
        return len(self.kids) + (1 if type(self.parent) is Cell else 0)

    def __getitem__(self, key):
        return self.kids[key]
    def __setitem__(self, key, value):
        self.kids[key] = value

    def get_vertices(self, index):
        vertices = []
        new_index = index + 1
        for kid in self.kids:
            vertices.append('{i} {j} 1'.format(i=index, j=new_index))
            vertices += kid.get_vertices(new_index)

            new_index += kid.total_cells()

        return vertices

    def export(self, filename):
        vertices = self.get_vertices(1)

        with open(filename, 'w') as f:
            f.write('*Vertices {n}\n\n'.format(n=len(vertices) + 1))
            f.write('*Arcs\n')
            for vertex in vertices:
                f.write(vertex + '\n')

    def reroot(self, base_nodes, parent=None):
        if parent is None:
            base_nodes.append(self)
        else:
            self.kids.remove(parent)

        if self.parent is None:
            base_nodes.remove(self)
        else:
            self.kids.append(self.parent)
            self.parent.reroot(base_nodes, self)

        self.parent = parent



def export(lst, filename):
    count = 1
    vertices = []

    for elem in lst:
        vertices += elem.get_vertices(count)
        count += elem.total_cells()

    with open(filename, 'w') as f:
        f.write('*Vertices {n}\n\n'.format(n=len(vertices) + len(lst)))
        f.write('*Arcs\n')
        for vertex in vertices:
            f.write(vertex + '\n')


max_size = 8*2048

class Population(object):
    """docstring for Population."""
    def __init__(self):
        global max_size
        self.pop = [[] for _ in range(max_size)]
        self.lm = 1

    def size(self, i=None):
        if i is None:
            return sum([len(self.pop[l]) for l in range(1, self.lm+1)])
        else:
            return len(self.pop[i])

    def total_cells(self):
        return sum([len(self.pop[l])*l for l in range(1, self.lm+1)])

    def lmax(self):
        return self.lm

    def update_lmax(self):
        global max_size
        lm_new = 0
        for l in range(1, self.lm+2):
            if len(self.pop[l])>0:
                lm_new = l
        self.lm = lm_new
        return self

    def add(self, cell):
        cell.parent = self
        size = cell.total_cells()
        self.pop[size].append(cell)
        return self

    def remove(self, cell, l=None):
        if l is None:
            l = cell.total_cells()
        try:
            self.pop[l].remove(cell)
        except ValueError:
            print("cellID: ", cell)
            print("current size: ", cell.total_cells())
            print("previous size: ", l)
            print("\nelements list (i, size, list)")
            for i in range(10):
                print(i, len(self.pop[i]), self.pop[i])

            print()
            traceback.print_exc(file=sys.stdout)
            exit()
        return self

    def pop_random(self, l):
        element = randrange(len(self.pop[l]))
        return self.pop[l].pop(element)

    def transfer(self, cell, old_l):
        self.remove(cell, old_l)
        self.add(cell)

    def __getitem__(self, key):
        return self.pop[key]
    def __setitem__(self, key, value):
        self.pop[key] = value


REP = 0
DEA = 1

def R1(C, p):
    x = 1.0
    for _ in range(10):
        x = x - (exp((1.-p)*x*x)*p*p + 4.*C*x*(-3. + 2.*(1.-p)*x*x))/(-2.*((p-1.)*x*exp((1.-p)*x*x)*p*p + 6.*C*(1. + 2.*(p-1.)*x*x)))
    return x*x*(1. - 8.*C/(p*p)*x*exp((p-1.)*x*x))


class Rates(object):
    """docstring for Rates."""
    def __init__(self, k_reproduction, k_death, stickiness, C, p):
        global max_size
        self.k_reproduction = k_reproduction
        self.k_death        = k_death
        self.stickiness     = stickiness

        self.reproduction = [0.0 for _ in range(max_size)]
        self.death        = [0.0 for _ in range(max_size)]

        self.fitness      = [0.0 for _ in range(max_size)]
        self.fitness[1]   = R1(C, p)
        for i in range(2, p):
            self.fitness[i] = i*self.fitness[1]
        for i in range(p, max_size):
            k = i%p
            m = i//p
            self.fitness[i] = i*pow(p,4.)*pow(m+1.,4.*k/p)*pow(m,4.*(p-k)/p)/(pow(i,4.)*432.*C*C)

        self.total_rate = 0.0

    def update(self, population):
        global max_size
        self.total_rate = 0.0

        lmax = population.lmax()
        for l in range(1, lmax+1):
            rate = self.k_reproduction*self.fitness[l]*population.size(l)
            self.reproduction[l] = rate
            self.total_rate += rate

        total_cells = population.total_cells()
        for l in range(1, lmax+1):
            rate = self.k_death*l*population.size(l)*total_cells
            self.death[l] = rate
            self.total_rate += rate

    def choose_reaction(self, lmax):
        global REP, DEA

        threshold = random()*self.total_rate
        rate = 0.0

        for l in range(1, lmax+1):
            rate += self.reproduction[l]
            if rate > threshold:
                return (REP, l)

        for l in range(1, lmax+1):
            rate += self.death[l]
            if rate > threshold:
                return (DEA, l)

        print("Error: could not find a reaction")
        print(self.total_rate, threshold)

def update(population, next_reaction, l, split_probability):
    global REP, DEA, k

    if next_reaction == REP:
        if random() < split_probability:
            population.add(Cell(k, population))
        else:
            group = population.pop_random(l)
            group.reproduce_random_member()
            population.add(group)
    elif DEA:
        group = population.pop_random(l)
        lst = group.kill_random_member()

        try:
            for group in lst:
                population.add(group)
        except TypeError:
            print("number of neighbours: ", len(group))
            traceback.print_exc(file=sys.stdout)
            exit()
    else:
        print("Unrecognized reaction type:", next_reaction)
        exit()

    population.update_lmax()

def save(t):
    return not (t%1000)

def save_results(f, t, population, reaction, l):
    lav  = population.total_cells()/population.size()
    lmax = population.lmax()

    f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(t, 0, lmax, lav, population.size(), population.total_cells(), reaction, l))


split_probability = 1.0-stickiness*stickiness


rates = Rates(k_rep, k_dea, stickiness, C, p)

population = Population()
for _ in range(initial_population):
    population.add(Cell(k, population))


from datetime import datetime
import socket

filename = "results/data_z{}_p{}_kr{}_kd{}_st{}_c{}.dat".format(k, p, k_rep, k_dea, stickiness, C)

with open(filename, "w") as f:
    f.write("##################################################################\n"   )
    f.write("# General data file for snowflake group with z = {}\n".format(k)        )
    f.write("#\n"                                                                    )
    f.write("# Run at {}\n".format(str(datetime.now()))                              )
    f.write("#     on {}\n".format(socket.gethostname())                             )
    f.write("#\n"                                                                    )
    f.write("# Parameters:\n"                                                        )
    f.write("#     k_reprod    {}\n".format(k_rep)                                   )
    f.write("#     k_death     {}\n".format(k_dea)                                   )
    f.write("#     stickiness  {}\n".format(stickiness)                              )
    f.write("#     p           {}\n".format(p)                                       )
    f.write("#     c           {}\n".format(C)                                       )
    f.write("##################################################################\n\n" )
    f.write("# iterate\ttime\tlmax\taverage_size\tgroup_number\ttotal_cells_number\treaction\tl\n")
    
    for t in range(t_max):
        rates.update(population)
        next_reaction, l = rates.choose_reaction(population.lmax())
        update(population, next_reaction, l, split_probability)

        if save(t):
            save_results(f, t, population, next_reaction, l)

