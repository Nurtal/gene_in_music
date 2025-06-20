from pydub import AudioSegment


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
    


if __name__ == "__main__":

    audio_file_list = ['signals/aorta/GTEX-1HSMP-0826-SM-A9SKW_signal.wav', 'signals/aorta/GTEX-ZF3C-1426-SM-4WWCD_signal.wav']


    concatenate_audio_file(audio_file_list, "/tmp/zog.wav")
    
