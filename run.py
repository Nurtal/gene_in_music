# usual import
import os
import glob

# import module
import craft_toy_data
import build_signal
import simple_clf
import craft_data


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


def simple_run():
    """Simple binary classification on tissue dataset, use random gene order and random
    gene selection, basically there just to make sure this stuff compile on real data"""

    # params
    n_genes = 100
    audio_diration = 4.0

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
    fil_list_a = glob.glob("signals/aorta/*.wav")
    fil_list_a = glob.glob("signals/coronary/*.wav")

    # run classification
    simple_clf.run_svm_clf(file_list_a, file_list_b)

if __name__ == "__main__":

    # toy_run()
    simple_run()
