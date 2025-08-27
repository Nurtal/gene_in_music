import glob
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, adjusted_rand_score, normalized_mutual_info_score
import pandas as pd



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
    
    


def evaluate_clustering(prediction_file:str, manifest_file:str, result_file:str) -> None:
    """Run clustering evaluation
    
    Args:
        - prediction_file (str) : csv file containing cluster prediction
        - manifest_file (str) : csv file containg true label of files
        - result_file (str) : output file to save the evaluation metrics
    
    """

    # load file to prediction
    file_to_prediction = {}
    df = pd.read_csv(prediction_file)
    for index, row in df.iterrows():
        file_to_prediction[row['FILE']] = row['LABEL']

    # load file to label    
    file_to_label = {}
    df = pd.read_csv(manifest_file)
    for index, row in df.iterrows():
        file_to_label[row['FILE']] = row['LABEL']

    # extract pred and true
    y_true = []
    y_pred = []
    for k in list(file_to_label.keys()):
        y_true.append(file_to_label[k])
        y_pred.append(file_to_prediction[k])

    # Compute adjusted Rand Index
    ari = adjusted_rand_score(y_true, y_pred)

    # Compute normalized Mutual Information
    nmi = normalized_mutual_info_score(y_true, y_pred)

    # save results
    output_data = open(result_file, "w")
    output_data.write("METRIC,VALUE\n")
    output_data.write(f"ADJUSTED-RAND-INDEX,{ari}\n")
    output_data.write(f"NORMALIZED-MUTUAL-INFORMATION,{nmi}\n")
    output_data.close()




if __name__ == "__main__":


    # prepare data & manifest
    file_list = []
    file_list_a = glob.glob("demo/group_a/*.wav")
    file_list_b = glob.glob("demo/group_b/*.wav")
    manifest_file = open("demo/manifest.csv", "w")
    manifest_file.write("FILE,LABEL\n")
    for elt in file_list_a:
        file_list.append(elt)
        manifest_file.write(f"{elt},0\n")
    for elt in file_list_b:
        file_list.append(elt)
        manifest_file.write(f"{elt},1\n")
    manifest_file.close()
    

    # run_kmeans(file_list, 2, 4, "/tmp/results.csv", 4.0, 2)

    evaluate_clustering("/tmp/results.csv", "demo/manifest.csv", "/tmp/cluster.csv")    
