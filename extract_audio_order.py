from pydub import AudioSegment
import glob


def concatenate_audio_file(audio_file_list:list, audio_file_output:str) -> None:
    """Concatenate audio file in the order presented in audio_file_list, save concatenate track to audio_file_output
    
    Args:
        audio_file_list (list) : list of .wav files
        audio_file_output (str) : path to output file
    
    """

    audio_concatene = AudioSegment.from_wav(audio_file_list[0])
    for f in audio_file_list[1:]:
        audio = AudioSegment.from_wav(f)
        audio_concatene += audio
    
    audio_concatene.export(audio_file_output, format="wav")


def extract_manifest(signals_folder:str) -> dict:
    """Extract id for which we have all pathways, return a dictionnary with class as key and associated ids as values

    Args:
        signals_folder (str) : path to signal folder, should be structures as follows : pathway/class/wav_file

    Returns:
        (dict) : class to list of associated ids for which we have a wav file for each pathway
    
    """

    # screen signals folder
    pathway_to_class_to_id = {}
    class_list = []
    pathway_list = []
    for pathway_folder in glob.glob(f"{signals_folder}/*"):
        pathway_name = pathway_folder.split("/")[-1]
        pathway_to_class_to_id[pathway_name] = {}
        if pathway_name not in pathway_list:
            pathway_list.append(pathway_name)
        for class_folder in glob.glob(f"{pathway_folder}/*"):
            class_name = class_folder.split("/")[-1]
            pathway_to_class_to_id[pathway_name][class_name] = []
            if class_name not in class_list:
                class_list.append(class_name)
            for signal_file in glob.glob(f"{class_folder}/*.wav"):
                id_file = signal_file.split("/")[-1].replace("_signal.wav", "")
                pathway_to_class_to_id[pathway_name][class_name].append(id_file)

    # check presence of id in all pathways for each class
    class_to_file_to_keep = {}
    for class_name in class_list:
        class_to_file_to_keep[class_name] = []
        original_list = pathway_to_class_to_id[pathway_list[0]][class_name]
        for pathway_name in pathway_list[1:]:
            for i in pathway_to_class_to_id[pathway_name][class_name]:
                if i in original_list and i not in class_to_file_to_keep[class_name]:
                    class_to_file_to_keep[class_name].append(i)

    return class_to_file_to_keep

    

if __name__ == "__main__":

    audio_file_list = ['signals/aorta/GTEX-1HSMP-0826-SM-A9SKW_signal.wav', 'signals/aorta/GTEX-ZF3C-1426-SM-4WWCD_signal.wav']

    # concatenate_audio_file(audio_file_list, "/tmp/zog.wav")
    m = extract_manifest("data/signals_small")
    print(m)
    
