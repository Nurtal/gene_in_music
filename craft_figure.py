import pandas as pd
import matplotlib.pyplot as plt

def plot_signal(signal_file:str, image_file:str) -> None:
    """Plot signal from a signal file

    Args:
        - signal_file (str) : path to the signal_file
        - image_file (str) : path to the figure file to generate
    
    """

    # load data
    df = pd.read_csv(signal_file)
    x = list(df['x'])
    y = list(df['y'])

    # craft plot
    plt.figure(figsize=(8, 4))
    plt.bar(x, y, width=1, color='blue', alpha=0.8)

    # Drop legends and stuff
    plt.xticks([])  
    plt.yticks([])
    plt.gca().spines['top'].set_visible(False)  
    plt.gca().spines['right'].set_visible(False)  
    plt.gca().spines['left'].set_visible(False)  
    plt.gca().spines['bottom'].set_visible(False)  
    plt.axis("off")

    # save image
    plt.savefig(image_file)




if __name__ == "__main__":

    plot_signal("signals/136_signal.csv", "test.png")
