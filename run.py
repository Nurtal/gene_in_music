# usual import
import os
import glob

# import module
import craft_toy_data
import build_signal
import simple_clf


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
    build_signal.build_random_signal("data/toy_data.csv")
    
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



if __name__ == "__main__":

    toy_run()
