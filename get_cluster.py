import glob
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, adjusted_rand_score, normalized_mutual_info_score



# import local module
import extract_features


def run_kmeans(file_list:list, J:int, Q:int, result_file:str, audio_duration:float, k:int) -> None:
    """Run kmeans clustering on audio file
    
    Args:
        - file_list (list) : list of audio file to treat
        - J (int) : scat features parameters 1
        - Q (int) : scat features parameters 2
        - result_save (str) : path to the file for saving results
        - audio_duration (float) : duration of the audio samples (seconds)
        - k (int) : number of cluster to hunt
    
    """

    # load data
    X = []
    for fl in file_list:
        x = extract_features.extract_features(fl, J, Q)
        features = x.numpy().flatten()
        X.append(features)

    # run kmeans
    kmeans = KMeans(n_clusters=k, random_state=42, n_init="auto")
    y_pred = kmeans.fit_predict(X)

    # save prediction
    i = 0
    prediction = open(result_file, "w")
    prediction.write("FILE,LABEL\n")
    for f in file_list:
        prediction.write(f"{f},{y_pred[i]}\n")
        i+=1
    prediction.close()

    # compute & save silhouette score
    sil_score = silhouette_score(X, y_pred)
    metric_file = open(f"{result_file.replace('.csv', '.log')}", "w")
    metric_file.write(f"Silouhette_score : {sil_score}\n")
    metric_file.close()
    
    


def evaluate_clustering(prediction_file, manifest_file):
    """ """


    sil_score = silhouette_score(X, y_pred)



if __name__ == "__main__":


    file_list = glob.glob("demo/group_a/*.wav")
    for elt in glob.glob("demo/group_b/*.wav"):
        file_list.append(elt)

    run_kmeans(file_list, 2, 4, "/tmp/results.csv", 4.0, 2)

    
