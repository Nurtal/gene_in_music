import pandas as pd
import random
import shutil
import os
import glob
import yaml

# local importation
import build_signal
import extract_features
import simple_clf


def generate_individual(genes:list, possible_positions:list) -> dict:
    """Create one random individual

    Args:
        - genes (list) : list of genes extracted from data file
        - possible_poistions (list) : list of possible positions to adopt

    Returns:
        - (dict) : gene to position
    
    """
    # randomize position order
    random.shuffle(possible_positions)

    # craft and return individual
    return {gene: pos for gene, pos in zip(genes, possible_positions)}

def generate_population(size:int, genes:list, possible_positions:list) -> list:
    """Generate a list of individual

    Args:
        - size (int) : n individual in pop
        - genes (list) : list of genes extracted from data file
        - possible_poistions (list) : list of possible positions to adopt

    Returns:
        - (list) : list of individuals, i.e a population
    
    """
    return [generate_individual(genes, possible_positions) for _ in range(size)]



def fitness(gene_to_pos, result_folder, configuration_file):
    """compute fitness for a given individual (gene_to_pos), in this case use auc from log classification"""

    # load configuration
    with open(configuration_file, "r") as f:
        config = yaml.safe_load(f)

    # load data
    df = pd.read_csv(config['data_file'])

    #--------------#
    # Prepare Data #
    #--------------#

    # clean result folder
    if os.path.isdir(result_folder):
        shutil.rmtree(result_folder)
    os.mkdir(result_folder)
        
    # cleaning data
    var_to_keep = ['ID']
    label_list = []
    for elt in gene_to_pos.keys():
        var_to_keep.append(elt)
    for label in df['GROUP']:
        if label not in label_list:
            label_list.append(label)
    for label in label_list:
        df_grp = df[df['GROUP'] == label]
        df_grp = df_grp[var_to_keep]
        df_grp.to_csv(f"{result_folder}/data_group_{label}.csv", index=False)

    # build signal
    for label in label_list:
        build_signal.build_signal_from_computed_positions(
                                                          f"{result_folder}/data_group_{label}.csv",
                                                          f"{result_folder}/group_{label}",
                                                          gene_to_pos
                                                      )

    # turn into audio files
    for label in label_list:
        for signal_file in glob.glob(f"{result_folder}/group_{label}/*.csv"):
            build_signal.turn_signal_into_audio(signal_file, config['audio_duration'])
    
    # prepare data for classification
    label_to_audio_list = {}
    for label in label_list:
        label_to_audio_list[label] = glob.glob(f"{result_folder}/group_{label}/*.wav")

    # take samples
    for label in label_list:
        audio_file = label_to_audio_list[label][random.randint(0, len(label_to_audio_list[label])-1)]
        extract_features.display_features(audio_file, config['J'], config['Q'], f"{result_folder}/signal_sample_group_{label}.png")        

    # ---------------#
    # Run Classifier #
    #----------------#

    # deal with binary log
    auc = 0
    if config['classifier'] == 'log':
        if len(label_list) == 2:
            auc = simple_clf.run_log_clf(
                    label_to_audio_list[label_list[0]],
                    label_to_audio_list[label_list[1]],
                    config['J'],
                    config['Q'],
                    f"{result_folder}/results.csv",
                    config['audio_duration']
            )
        else:
            print("[!] Can't run binary claffication with n labels != 2")

    return auc


def selection(population:list, result_folder:str, configuration_file:str) -> dict:
    """Sélection (tournoi à 2)
    For some reason classifier sometime crash because of nan in data, if its the case fitness score is
    condisedered 0
    
    """
    
    # pick random candidates
    a, b = random.sample(population, 2)

    # compute score a
    score_a = 0
    try:
        score_a = fitness(a, result_folder, configuration_file)
    except:
        pass

    # compute score a
    score_b = 0
    try:
        score_b = fitness(b, result_folder, configuration_file)
    except:
        pass
    
    return (a, score_a) if score_a > score_b else (b, score_b)


def evaluate_population(population:list, result_folder:str, configuration_file:str):
    """ """

    id_to_score = {}
    id_to_ind = {}
    id = 0
    for ind in population:
        score = 0
        id +=1
        try:
            score = fitness(ind, result_folder, configuration_file)
        except:
            pass

        id_to_score[id] = score
        id_to_ind[id] = ind

    return (id_to_score, id_to_ind)
        


def crossover(parent1, parent2, genes):
    """ """

    # On extrait juste les permutations
    p1_perm = [parent1[g] for g in genes]
    p2_perm = [parent2[g] for g in genes]

    size = len(p1_perm)
    start, end = sorted(random.sample(range(size), 2))

    # on prend un segment du parent1
    child_perm = [None] * size
    child_perm[start:end] = p1_perm[start:end]

    # on complète avec l'ordre du parent2
    fill_values = [x for x in p2_perm if x not in child_perm]
    idx = 0
    for i in range(size):
        if child_perm[i] is None:
            child_perm[i] = fill_values[idx]
            idx += 1

    return {gene: pos for gene, pos in zip(genes, child_perm)}




def mutate(individual, mutation_rate, genes):
    if random.random() < mutation_rate:
        g1, g2 = random.sample(genes, 2)
        individual[g1], individual[g2] = individual[g2], individual[g1]









def extract_optimal_position(data_file, result_file):
    """ """

    # load data & set paramaters
    df = pd.read_csv(data_file)
    genes = list(df.keys())[1:-1]
    possible_positions = list(range(len(genes)))  # autant de positions que de gènes
    pop_size = 2
    n_generation = 1
    mutation_rate = 0.3
    result_folder = "/tmp/ga_gim"
    configuration_file = "ressources/example_config.yaml"
    result_data = []

    # init population
    perfect_order = {}
    best_score = 0
    best_pos = None
    population = generate_population(pop_size, genes, possible_positions)
    for gen in range(n_generation):
        new_population = []
        for _ in range(pop_size):

            # take random parents
            s1 = selection(population, result_folder, configuration_file)
            parent1 = s1[0]
            score1 = s1[1]
            s2 = selection(population, result_folder, configuration_file)
            parent2 = s2[0]
            score2 = s2[1]

            # look for perfect parent
            if score1 == 1.0:
                print("FIND PERFECT ORDER")
                best_pos = parent1
                best_score = score1
                break
            if score2 == 1.0:
                print("FIND PERFECT ORDER")
                best_pos = parent2
                best_score = score2
                break

            # create child
            child = crossover(parent1, parent2, genes)
            mutate(child, mutation_rate, genes)
            new_population.append(child)
        
        # assemble new population
        population = new_population
        results = evaluate_population(population, result_folder, configuration_file)
        id_to_score = results[0]
        id_to_pos = results[1] 

        # update data results
        for id in id_to_score:
            result_data.append(
                {
                    "GENERATION": gen,
                    "POSITION":id_to_pos[id],
                    "SCORE" : id_to_score[id]
                }
            )
            
            # find best positions
            if id_to_score[id] > best_score:
                best_score = id_to_score[id]
                best_pos = id_to_pos[id]

        # check if best position is perfect positions
        if best_score == 1.0:
            print("FIND PERFECT ORDER")
            break

    # save data results
    df = pd.DataFrame.from_dict(result_data)
    df.to_csv(result_file, index=False)

    # return best results
    return (best_pos, best_score)


if __name__ == "__main__":


    m = extract_optimal_position("data/fake_gene_data.csv", "/tmp/ga_results.csv")
    print(m)
    
