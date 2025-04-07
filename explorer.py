from itertools import product
from tqdm import tqdm


# import module
import run



def run_binary_exploration():
    """Explore different combination"""

    # params
    param_J = [1, 2,3,4,5,6,7,8,9,10]
    param_Q = [1, 2,3,4,5,6,7,8,9,10]
    param_duration = [1.0, 2.0, 3.0,4.0,5.0,6.0,7.0]

    # Get all possible combination
    combinations = list(product(param_J, param_Q, param_duration))
    
    # prepare result folder
    if os.path.exists("exploration") and os.path.isdir("exploration"):
        shutil.rmtree("exploration")
    os.mkdir("exploration")

    # Affiche ou utilise
    cmpt = 0
    for combo in tqdm(combinations, desc="Test des combinaisons"):

        # extract params
        J = combo[0]
        Q = combo[1]
        audio_duration = combo[2]
        output_folder = f"exploration/try_{cmpt}"
        result_file = f"exploration/try_{cmpt}/results.csv"
        preprocess_data = True

        if cmpt > 0:
            preprocess_data = False

        # run
        run.simple_binary_run(output_folder,
                              preprocess_data,
                              audio_duration,
                              J,
                              Q,
                              result_file
                          )

        # update cmpt
        cmpt +=1





if __name__ == "__main__":
    
    run_binary_exploration()
