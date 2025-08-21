import random
import pandas as pd
import os


def craft_toy_data(nb_patient_group_a:int, nb_patient_group_b:int, nb_noisy_genes:int) -> None:
    """Craft a toy dataset with 2 groups describe by 10 'genes'.
    Save result file in data subfolder, create it of not exist

    Args:
        - nb_patient_group_a (int) : number of patient in group a
        - nb_patient_group_b (int) : number of patient in group b
        - nb_noisy_genes (int) : number of noisy genes to add
    
    """

    # vars
    data = []
    cmpt = 0

    # craft vectors - group a
    for x in range(nb_patient_group_a):

        # forge vector
        vector = {
            "ID":cmpt,
            "A1":random.randint(10,20),
            "A2":random.randint(10,20),
            "A3":random.randint(10,20),
            "B1":random.randint(0,10),
            "B2":random.randint(0,10),
            "B3":random.randint(0,10),
            "C1":random.randint(0,20),
            "C2":random.randint(0,20),
            "C3":random.randint(0,20),
            "D1":random.randint(0,100),
        }

        # add noisy gene
        for x in range(nb_noisy_genes):
            vector[f"N{x}"] = random.randint(8,12)

        # update cmpt
        cmpt +=1

        # add vector to data
        data.append(vector)
    
    # craft vectors - group b
    for x in range(nb_patient_group_b):

        # forge vector
        vector = {
            "ID":cmpt,
            "A1":random.randint(0,10),
            "A2":random.randint(0,10),
            "A3":random.randint(0,10),
            "B1":random.randint(10,20),
            "B2":random.randint(10,20),
            "B3":random.randint(10,20),
            "C1":random.randint(0,20),
            "C2":random.randint(0,20),
            "C3":random.randint(0,20),
            "D1":random.randint(0,100),
        }

        # add noisy gene
        for x in range(nb_noisy_genes):
            vector[f"N{x}"] = random.randint(8,12)

        # update cmpt
        cmpt +=1

        # add vector to data
        data.append(vector)
    

    # craft dataframe
    df = pd.DataFrame(data)

    # save dataframe
    if not os.path.isdir("data"):
        os.mkdir("data")
    df.to_csv("data/toy_data.csv", index=False)



def craft_fake_gene_dataset(nb_patient_group_a:int, nb_patient_group_b:int) -> None:
    """Craft a toy dataset with 2 groups describe by genes from 2 distinct pathways : IFN et mTOR.
    Save result file in data subfolder, create it of not exist

    Args:
        - nb_patient_group_a (int) : number of patient in group a
        - nb_patient_group_b (int) : number of patient in group b
    """

    # gene sets
    cmpt = 0
    data = []
    ifn_list = [
        "IFNA1",
        "IFNA2",
        "IFNB1",
        "STAT1",
        "STAT2",
        "IRF7",
        "IRF9",
        "ISG15",
        "MX1",
        "OAS1",
        "IFIT1",
        "IFNG",
        "JAK1",
        "JAK2",
        "CXCL10",
        "GBP1"
    ]
    mtor_list = [
        "MTOR",
        "RPTOR",
        "RICTOR",
        "AKT1",
        "TSC1",
        "TSC2",
        "S6K1",
        "4EBP1"
    ]

    # craft grouop a
    for x in range(nb_patient_group_a):

        # update cmpt
        cmpt +=1

        # init vector
        vector = {"ID":cmpt}

        # add ifn
        for g in ifn_list:
            vector[g] = random.randint(5,10)

        # add mtor
        for g in mtor_list:
            vector[g] = random.randint(0,10)

        # add vector to data
        data.append(vector)
    
    
    # craft grouop b
    for x in range(nb_patient_group_b):

        # update cmpt
        cmpt +=1

        # init vector
        vector = {"ID":cmpt}

        # add ifn
        for g in ifn_list:
            vector[g] = random.randint(0,5)

        # add mtor
        for g in mtor_list:
            vector[g] = random.randint(0,10)

        # add vector to data
        data.append(vector)


    # craft dataframe
    df = pd.DataFrame(data)

    # save dataframe
    if not os.path.isdir("data"):
        os.mkdir("data")
    df.to_csv("data/fake_gene_data.csv", index=False)


if __name__ == "__main__":

    craft_toy_data(100, 100, 25)
    craft_fake_gene_dataset(50,50)
