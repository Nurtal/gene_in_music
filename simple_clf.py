import glob
import extract_features
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
import pandas as pd
import umap.umap_ as umap


def run_svm_clf(file_list_1:list, file_list_2:list):
    """
    Simple exemple case, extract features from scats and used it to train a SVM
    """

    # params
    J = 6 
    Q = 8
    duration = 16000

    X = []
    y = []
    for fl in file_list_1:
        x = extract_features.extract_features(fl, J, Q)
        features = x.numpy().flatten()
        X.append(features)
        y.append("class_a")
    for fl in file_list_2:
        x = extract_features.extract_features(fl, J, Q)
        features = x.numpy().flatten()
        X.append(features)
        y.append("class_b")
    
    # Division en ensemble d'entraînement et de test
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # Entraînement du modèle SVM
    clf = SVC(kernel="linear", C=1.0)
    clf.fit(X_train, y_train)

    # Prédiction et évaluation
    y_pred = clf.predict(X_test)
    accuracy = accuracy_score(y_test, y_pred)
    print(f"[CLF][SVM] ACC : {accuracy * 100:.2f}%")
    
  

def run_log_clf(file_list_1:list, file_list_2:list, J:int, Q:int, result_file:str, audio_duration:float) -> float:
    """
    Simple exemple case, extract features from scats and used it to train a logistic regression

    Args:
        - file_list_1 (list) : list of file for the class a
        - file_list_2 (list) : list of file for the class b
        - J (int) : scat features parameters 1
        - Q (int) : scat features parameters 2
        - result_save (str) : path to the file for saving results
        - audio_duration (float) : duration of the audio samples (seconds)

    Returns:
        - (float) : auc
    
    """

    # load data
    X = []
    y = []
    for fl in file_list_1:
        x = extract_features.extract_features(fl, J, Q)
        features = x.numpy().flatten()
        X.append(features)
        y.append("class_a")
    for fl in file_list_2:
        x = extract_features.extract_features(fl, J, Q)
        features = x.numpy().flatten()
        X.append(features)
        y.append("class_b")
    
    # Division en ensemble d'entraînement et de test
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # Entraînement du modèle SVM
    clf = LogisticRegression()
    clf.fit(X_train, y_train)

    # Prédiction et évaluation
    y_pred = clf.predict(X_test)
    accuracy = accuracy_score(y_test, y_pred)
    print(f"[CLF][LOG-REF] ACC : {accuracy * 100:.2f}%")

    # compute AUC
    y_probs = clf.predict_proba(X_test)[:, 1]
    auc = roc_auc_score(y_test, y_probs)
    print(f"[CLF][LOG-REF] AUC : {auc}")

    # save results
    output_file = open(result_file, "w")
    output_file.write("METRIC,VALUE\n")
    output_file.write("CLF,Logistic-Regression\n")
    output_file.write(f"J,{J}\n")
    output_file.write(f"Q,{Q}\n")
    output_file.write(f"Audio-Duration,{audio_duration}\n")
    output_file.write(f"ACC,{accuracy}\n")
    output_file.write(f"AUC,{auc}\n")
    output_file.close()

    # return auc
    return auc



def run_direct_log_clf(data_file_a:str, data_file_b:str, result_file:str) -> None:
    """
    Simple exemple case, train a logistic regression directly on the rnaseq data

    Args:
        - data_file_a (str) : data file containing observations for class a (csv file)
        - data_file_b (str) : data file containing observations for class b (csv file)
        - result_save (str) : path to the file for saving results
    
    """

    # load data
    X = []
    y = []
    df_a = pd.read_csv(data_file_a)
    df_a = df_a.drop(columns=['ID'])
    for index, row in df_a.iterrows():
        vector = []
        for k in row.keys():
            vector.append(float(row[k]))
        X.append(vector)
        y.append("class_a")
    df_b = pd.read_csv(data_file_b)
    df_b = df_b.drop(columns=['ID'])
    for index, row in df_b.iterrows():
        vector = []
        for k in row.keys():
            vector.append(float(row[k]))
        X.append(vector)
        y.append("class_b")

    # Division en ensemble d'entraînement et de test
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # Entraînement du modèle SVM
    clf = LogisticRegression()
    clf.fit(X_train, y_train)

    # Prédiction et évaluation
    y_pred = clf.predict(X_test)
    accuracy = accuracy_score(y_test, y_pred)
    print(f"[CLF][LOG-REF] ACC : {accuracy * 100:.2f}%")

    # compute AUC
    y_probs = clf.predict_proba(X_test)[:, 1]
    auc = roc_auc_score(y_test, y_probs)
    print(f"[CLF][LOG-REF] AUC : {auc}")

    # save results
    output_file = open(result_file, "w")
    output_file.write("METRIC,VALUE\n")
    output_file.write("CLF,Logistic-Regression\n")
    output_file.write(f"ACC,{accuracy}\n")
    output_file.write(f"AUC,{auc}\n")
    output_file.close()



def run_umap_log_clf(data_file_a:str, data_file_b:str, result_file:str) -> None:
    """
    Train a logistic regression on a umap crafted from rnaseq data

    Args:
        - data_file_a (str) : data file containing observations for class a (csv file)
        - data_file_b (str) : data file containing observations for class b (csv file)
        - result_save (str) : path to the file for saving results
    
    """

    # load data
    X = []
    y = []
    df_a = pd.read_csv(data_file_a)
    df_a = df_a.drop(columns=['ID'])
    for index, row in df_a.iterrows():
        vector = []
        for k in row.keys():
            vector.append(float(row[k]))
        X.append(vector)
        y.append("class_a")
    df_b = pd.read_csv(data_file_b)
    df_b = df_b.drop(columns=['ID'])
    for index, row in df_b.iterrows():
        vector = []
        for k in row.keys():
            vector.append(float(row[k]))
        X.append(vector)
        y.append("class_b")

    # Convertir en numpy array pour UMAP
    data_array = np.array(X)

    # Initialiser UMAP (n_components=2 pour une projection 2D)
    reducer = umap.UMAP(n_components=2, random_state=42)

    # Appliquer la réduction de dimension
    X = reducer.fit_transform(data_array)

    # Division en ensemble d'entraînement et de test
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # Entraînement du modèle SVM
    clf = LogisticRegression()
    clf.fit(X_train, y_train)

    # Prédiction et évaluation
    y_pred = clf.predict(X_test)
    accuracy = accuracy_score(y_test, y_pred)
    print(f"[CLF][LOG-REF] ACC : {accuracy * 100:.2f}%")

    # compute AUC
    y_probs = clf.predict_proba(X_test)[:, 1]
    auc = roc_auc_score(y_test, y_probs)
    print(f"[CLF][LOG-REF] AUC : {auc}")

    # save results
    output_file = open(result_file, "w")
    output_file.write("METRIC,VALUE\n")
    output_file.write("CLF,Logistic-Regression\n")
    output_file.write(f"ACC,{accuracy}\n")
    output_file.write(f"AUC,{auc}\n")
    output_file.close()



def run_log_binary_clf_on_audio():
    """ """


if __name__ == "__main__":


    # run_svm_clf()

    # run_direct_log_clf("/tmp/zog/data/gene_reads_artery_aorta_HALLMARK_ADIPOGENESIS.csv", "/tmp/zog/data/gene_reads_artery_coronary_HALLMARK_ADIPOGENESIS.csv", "/tmp/machin.csv")
    run_umap_log_clf("/tmp/zog/data/gene_reads_artery_aorta_HALLMARK_ADIPOGENESIS.csv", "/tmp/zog/data/gene_reads_artery_coronary_HALLMARK_ADIPOGENESIS.csv", "/tmp/machin.csv")
