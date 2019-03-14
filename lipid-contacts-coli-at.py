import sys
sys.path.append("~/Documents/PART_II/USEFUL/scripts")
import MDAnalysis as mda
import numpy as np
from contacts import LipidContactGenerator
import argparse

parser = argparse.ArgumentParser('Lipid Contacts - selecting the lipid')

parser.add_argument(
	"-l",
	type = str,
	help = "lipid name"
)

parser.add_argument(
	"-f",
	type = str,
	help = "GRO file"
)

parser.add_argument(
	"-t",
	type = str,
	help = "XTC file"
)

parser.add_argument(
	"-o",
	type = str,
	help = "output name"
)

args = parser.parse_args()

U = mda.Universe(args.f,args.t)

generate = LipidContactGenerator(U)

contacts = generate.build_contacts(ligand_selection="resname " + args.l, frameskip=1, cutoff=3, KDTree=True)

contacts.aggregate(group_protein_by="resid",group_ligand_by="segid",aggregate_function=lambda x:x.max())

data = contacts.time_aggregate(aggregate_function=lambda x:np.sum(x.values())/contacts.n_frames)

data.to_dataframe().to_csv("{}_{}-contacts.csv".format(args.o, args.l))

contacts = generate.build_contacts(ligand_selection="resname " + args.l, frameskip=1, cutoff=3, KDTree=True)

contacts.aggregate(group_protein_by="resid",group_ligand_by="resname",aggregate_function=lambda x:x.max())

data = contacts.time_aggregate(aggregate_function=lambda x:np.sum(x.values())/contacts.n_frames)

data.to_dataframe().to_csv("{}_contacts.csv".format(args.o))

