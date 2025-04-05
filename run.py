import os
import glob

import craft_toy_data
import build_signal


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

    # run classification

    # display results




if __name__ == "__main__":

    toy_run()
