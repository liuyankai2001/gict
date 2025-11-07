import os
import sys
from path_search.balance import *
from path_search.data import Data
from path_search.pathways import *
from path_search.reactions import *
from path_search.enzymes import *
from path_search.rank import *
from path_search.networkexpansion import *
from config import project_config
from path_search.visualize import *

def get_path_search(target_compound):

    dat = Data(target_compound)
    dat.readParametersFile()

    # exclude unbalanced reactions if set to do so in the parameters
    if dat.exclude_unbalanced == 1:
        # calculate balance for each reaction from the compound formulas
        if dat.calculate_balance == 1:
            bal = Balance()
            # if not os.path.exists('../data/reaction_balance.csv'):
            if not os.path.exists(project_config.DATA_DIR / 'reaction_balance.csv'):
                bal.createBalanceFile()
        dat.getBalancedReactionsDf()

    # do network expansion if pathways_or_networkexp parameter = 'n'
    if dat.pathways_or_networkexp == 'n':
        net = Network(dat)  # load network
        nexp = NetworkExpansion(dat,
                                net.G)  # create instance of the network expansion class and pass the data and parameters
        nexp.runExpansion()  # run network expansion
        exit()

    # do pathway search if the network expansion was not requested
    if 1 in dat.stages:
        dat.findPrecursorAndTargetLCSBID()

        net = Network(dat)

        print("Size of the network loaded: ", len(net.G))

        path = Pathways(dat, net.G)
        path.findPathwaysSet()

    if 2 in dat.stages:
        rxn = Reactions(dat)
        rxn.getAllPairs()
        rxn.assignPairUidToPairs()
        rxn.writePathwaysAsReactionsPw()

    if 3 in dat.stages:
        enz = Enzymes(dat)
        enz.getAllReactions()
        enz.assignECtoReactions()
        enz.writePathwaysAsEnzymesPw()

    if 4 in dat.stages:
        rank = Rank(dat)
        rank.getPathwayRankForAll()

    if 5 in dat.stages:
        try:

            vis = Visualization(dat)
            vis.printCompoundImages()
            vis.drawPathways()
        except ImportError as e:
            print("Visualization could not be done due to the following error:")
            print(e)
            if e.msg == "No module named 'rdkit'":
                print("Please, install rdkit and use rdkit environment for running ARBRE with visualization")
                print("First install anaconda: https://docs.anaconda.com/anaconda/install/index.html")
                print("Then install rdkit as described here: https://www.rdkit.org/docs/Install.html")
            pass

    if os.path.exists(dat.pathway_file) \
            and os.path.exists(dat.pathway_reactions_file) \
            and os.path.exists(dat.pathway_enzymes_file) \
            and os.path.exists(dat.ranked_pathway_file):
        os.remove(dat.pathway_file)
        os.remove(dat.pathway_reactions_file)
        os.remove(dat.pathway_enzymes_file)