# usual import
import os
import glob
import shutil
import random

# import module
import craft_toy_data
import build_signal
import simple_clf
import craft_data
import extract_gene_order
import extract_features
import craft_report

def toy_run():
    """ """

    # clean signal folder
    old_files = glob.glob(os.path.join("signals", "*"))
    for f in old_files:
        if os.path.isfile(f):
            os.remove(f)

    # generate toy dataset
    craft_toy_data.craft_toy_data(50, 50, 25)    

    # build signal
    build_signal.build_random_signal("data/toy_data.csv", "signals")
    
    # turn into audio files
    for signal_file in glob.glob("signals/*.csv"):
        build_signal.turn_signal_into_audio(signal_file, 4.0)

    # prepare data for classification
    file_list_a = []
    file_list_b = []
    for audio_file in glob.glob("signals/*.wav"):
        nb = int(audio_file.split("/")[1].split("_")[0])
        if nb < 50:
            file_list_a.append(audio_file)
        else:
            file_list_b.append(audio_file)
            
    # run classification
    simple_clf.run_svm_clf(file_list_a, file_list_b)


def simple_random_run():
    """Simple binary classification on tissue dataset, use random gene order and random
    gene selection, basically there just to make sure this stuff compile on real data"""

    # params
    n_genes = 100
    audio_duration = 4.0

    # generate datasets from gcts
    craft_data.craft_reduce_datasets(["data/gene_reads_artery_aorta.gct", "data/gene_reads_artery_coronary.gct"], n_genes)

    # build signal
    build_signal.build_random_signal("data/gene_reads_artery_aorta.csv", "signals/aorta")
    build_signal.build_random_signal("data/gene_reads_artery_coronary.csv", "signals/coronary")
    
    # turn into audio files
    for signal_file in glob.glob("signals/aorta/*.csv"):
        build_signal.turn_signal_into_audio(signal_file, audio_duration)
    for signal_file in glob.glob("signals/coronary/*.csv"):
        build_signal.turn_signal_into_audio(signal_file, audio_duration)

    # prepare data for classification
    file_list_a = glob.glob("signals/aorta/*.wav")
    file_list_b = glob.glob("signals/coronary/*.wav")

    # run classification
    simple_clf.run_svm_clf(file_list_a, file_list_b)
    

def simple_reduced_run(output_folder):
    """Simple binary classification on tissue dataset, use gene order computed from correlation and random
    gene selection, basically there just to make sure this stuff compile on real data"""

    # params
    n_genes = 100
    audio_duration = 4.0
    n_random_pick = 3
    J = 2
    Q = 4

    # prepare result folder
    if os.path.exists(output_folder) and os.path.isdir(output_folder):
        shutil.rmtree(output_folder)
    os.mkdir(output_folder)
    os.mkdir(f"{output_folder}/signal_samples")

    # generate datasets from gcts
    craft_data.craft_reduce_datasets(["data/gene_reads_artery_aorta.gct", "data/gene_reads_artery_coronary.gct"], n_genes)

    # compute gene order
    extract_gene_order.get_proximity_from_data(['data/gene_reads_artery_aorta.csv', 'data/gene_reads_artery_coronary.csv'], "data/prox_matrix.csv")
    gene_to_pos = extract_gene_order.build_order_from_proximity("data/prox_matrix.csv")
    
    # build signal
    build_signal.build_signal_from_computed_positions("data/gene_reads_artery_aorta.csv", "signals/aorta", gene_to_pos)
    build_signal.build_signal_from_computed_positions("data/gene_reads_artery_coronary.csv", "signals/coronary", gene_to_pos)
    
    # turn into audio files
    for signal_file in glob.glob("signals/aorta/*.csv"):
        build_signal.turn_signal_into_audio(signal_file, audio_duration)
    for signal_file in glob.glob("signals/coronary/*.csv"):
        build_signal.turn_signal_into_audio(signal_file, audio_duration)

    # prepare data for classification
    file_list_a = glob.glob("signals/aorta/*.wav")
    file_list_b = glob.glob("signals/coronary/*.wav")

    # take a pick at random files from a
    random_pick_a = random.sample(file_list_a, n_random_pick)
    for audio_file in random_pick_a:
        save_file = audio_file.split("/")[-1].replace(".wav", "_class_a.png")
        extract_features.display_features(audio_file, J, Q, f"{output_folder}/signal_samples/{save_file}")        

    # take a pick at random files from b
    random_pick_b = random.sample(file_list_b, n_random_pick)
    for audio_file in random_pick_b:
        save_file = audio_file.split("/")[-1].replace(".wav", "_class_b.png")
        extract_features.display_features(audio_file, J, Q, f"{output_folder}/signal_samples/{save_file}")        

    # un classification
    simple_clf.run_log_clf(file_list_a, file_list_b, J, Q, f"{output_folder}/results.csv", audio_duration)



def simple_binary_run(output_folder:str, preprocess_data:bool, audio_duration:float, J:int, Q:int, result_file:str):
    """Simple binary classification on tissue dataset, use gene order computed from correlation, used for parameters exploration

    WARNING : too much memory usage
    
    
    Args:
        - output_folder (str) : path to the result folder
        - preprocess_data (bool) : if set to true, craft datasets and build signals, else, skip this part
        - audio_duration (float) : duration of the audio samples (seconds)
        - J (int) : scat features parameters 1
        - Q (int) : scat features parameters 2
        - result_save (str) : path to the file for saving results
    
    """

    # params
    n_random_pick = 3

    # prepare result folder
    if os.path.exists(output_folder) and os.path.isdir(output_folder):
        shutil.rmtree(output_folder)
    os.mkdir(output_folder)
    os.mkdir(f"{output_folder}/signal_samples")

    # run data preprocessing
    if preprocess_data:

        # generate datasets from gcts
        craft_data.craft_datasets(["data/gene_reads_artery_aorta.gct", "data/gene_reads_artery_coronary.gct"])

        # compute gene order
        extract_gene_order.get_proximity_from_data(['data/gene_reads_artery_aorta.csv', 'data/gene_reads_artery_coronary.csv'], "data/prox_matrix.csv")
        gene_to_pos = extract_gene_order.build_order_from_proximity("data/prox_matrix.csv")
    
        # build signal
        build_signal.build_signal_from_computed_positions("data/gene_reads_artery_aorta.csv", "signals/aorta", gene_to_pos)
        build_signal.build_signal_from_computed_positions("data/gene_reads_artery_coronary.csv", "signals/coronary", gene_to_pos)
    
    # turn into audio files
    for signal_file in glob.glob("signals/aorta/*.csv"):
        build_signal.turn_signal_into_audio(signal_file, audio_duration)
    for signal_file in glob.glob("signals/coronary/*.csv"):
        build_signal.turn_signal_into_audio(signal_file, audio_duration)

    # prepare data for classification
    file_list_a = glob.glob("signals/aorta/*.wav")
    file_list_b = glob.glob("signals/coronary/*.wav")

    # take a pick at random files from a
    random_pick_a = random.sample(file_list_a, n_random_pick)
    for audio_file in random_pick_a:
        save_file = audio_file.split("/")[-1].replace(".wav", "_class_a.png")
        extract_features.display_features(audio_file, J, Q, f"{output_folder}/signal_samples/{save_file}")        

    # take a pick at random files from b
    random_pick_b = random.sample(file_list_b, n_random_pick)
    for audio_file in random_pick_b:
        save_file = audio_file.split("/")[-1].replace(".wav", "_class_b.png")
        extract_features.display_features(audio_file, J, Q, f"{output_folder}/signal_samples/{save_file}")        

    # un classification
    simple_clf.run_log_clf(file_list_a, file_list_b, J, Q, result_file, audio_duration)



def simple_binary_gsea_run(output_folder:str, preprocess_data:bool, audio_duration:float, J:int, Q:int):
    """ Perform binary log classification on each of the dataset crafted with gsea analysis
    
    Args:
        - output_folder (str) : path to the result folder
        - preprocess_data (bool) : if set to true, craft datasets and build signals, else, skip this part
        - audio_duration (float) : duration of the audio samples (seconds)
        - J (int) : scat features parameters 1
        - Q (int) : scat features parameters 2
    
    """

    # params
    n_random_pick = 3

    # prepare result folder
    if os.path.exists(output_folder) and os.path.isdir(output_folder):
        shutil.rmtree(output_folder)
    os.mkdir(output_folder)
    os.mkdir(f"{output_folder}/signal_samples")
    os.mkdir(f"{output_folder}/data")
    os.mkdir(f"{output_folder}/signals")
    os.mkdir(f"{output_folder}/results")
    os.mkdir(f"{output_folder}/results_direct")
    os.mkdir(f"{output_folder}/results_umap")

    # run data preprocessing
    if preprocess_data:

        # generate datasets from gcts
        craft_data.craft_gsea_dataset(["data/gene_reads_artery_aorta.gct", "data/gene_reads_artery_coronary.gct"], "data/h.all.v2024.1.Hs.entrez.gmt", f"{output_folder}/data")

        # loop over generated gsea file
        for data_file in glob.glob(f"{output_folder}/data/gene_reads_artery_aorta*.csv")[5:6]: # TODO LOOP ON ALL LIST

            # extract gene set name & prepare output dirs
            gene_set = data_file.split("/")[-1].replace("gene_reads_artery_aorta_", "").replace(".csv", "")
            os.mkdir(f"{output_folder}/signals/{gene_set}")
            os.mkdir(f"{output_folder}/signals/{gene_set}/aorta")
            os.mkdir(f"{output_folder}/signals/{gene_set}/coronary")

            # get associated file
            associated_data_file = data_file.replace("gene_reads_artery_aorta", "gene_reads_artery_coronary")

            # compute gene order
            extract_gene_order.get_proximity_from_data([data_file, associated_data_file], f"{output_folder}/data/prox_matrix.csv")
            gene_to_pos = extract_gene_order.build_order_from_proximity(f"{output_folder}/data/prox_matrix.csv")
    
            # build signal
            build_signal.build_signal_from_computed_positions(data_file, f"{output_folder}/signals/{gene_set}/aorta", gene_to_pos)
            build_signal.build_signal_from_computed_positions(associated_data_file, f"{output_folder}/signals/{gene_set}/coronary", gene_to_pos)
    
            # turn into audio files
            for signal_file in glob.glob(f"{output_folder}/signals/{gene_set}/aorta/*.csv"):
                build_signal.turn_signal_into_audio(signal_file, audio_duration)
            for signal_file in glob.glob(f"{output_folder}/signals/{gene_set}/coronary/*.csv"):
                build_signal.turn_signal_into_audio(signal_file, audio_duration)

            # prepare data for classification
            file_list_a = glob.glob(f"{output_folder}/signals/{gene_set}/aorta/*.wav")
            file_list_b = glob.glob(f"{output_folder}/signals/{gene_set}/coronary/*.wav")

            # take a look at random files from a
            random_pick_a = random.sample(file_list_a, n_random_pick)
            for audio_file in random_pick_a:
                save_file = audio_file.split("/")[-1].replace(".wav", "_class_a.png")
                extract_features.display_features(audio_file, J, Q, f"{output_folder}/signal_samples/{save_file}")        

            # take a look at random files from b
            random_pick_b = random.sample(file_list_b, n_random_pick)
            for audio_file in random_pick_b:
                save_file = audio_file.split("/")[-1].replace(".wav", "_class_b.png")
                extract_features.display_features(audio_file, J, Q, f"{output_folder}/signal_samples/{save_file}")        

            # un classification
            result_file = f"{output_folder}/results/{gene_set}_log_clf.csv"
            simple_clf.run_log_clf(file_list_a, file_list_b, J, Q, result_file, audio_duration)

            # un classification - direct
            result_file = f"{output_folder}/results_direct/{gene_set}_log_clf.csv"
            simple_clf.run_direct_log_clf(data_file, associated_data_file, result_file)

            # un classification - umap
            result_file = f"{output_folder}/results_umap/{gene_set}_log_clf.csv"
            simple_clf.run_umap_log_clf(data_file, associated_data_file, result_file)

    # craft report
    craft_report.craft_run_report(output_folder)
        

if __name__ == "__main__":

    # toy_run()
    # simple_reduced_run("/tmp/zog")

    output_folder = "/tmp/zogzog"
    audio_duration = 10.0
    J = 2
    Q = 7
    
    simple_binary_gsea_run(output_folder, True, audio_duration, J, Q)
