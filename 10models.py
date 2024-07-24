import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.preprocessing import LabelEncoder, StandardScaler, RobustScaler
from sklearn.metrics import make_scorer, confusion_matrix
from sklearn.feature_selection import SelectKBest, f_classif, RFE
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.pipeline import Pipeline
from sklearn.utils import resample
from sklearn.decomposition import PCA
import warnings
import xgboost as xgb
import matplotlib.pyplot as plt
import seaborn as sns

# Suppress warnings in the output
warnings.filterwarnings('ignore')

# Load data from a CSV file
data = pd.read_csv('expr_inter_ad.csv')
X = data.drop('Group', axis=1)  # Features
y = data['Group'].replace({'CN': 0, 'Dementia': 1})  # Labels

# Encode labels to integers
le = LabelEncoder()
y = le.fit_transform(y)

# Determine majority and minority classes
majority_class = np.argmax(np.bincount(y))
minority_class = 1 - majority_class

# Get indices for majority and minority classes
majority_indices = np.where(y == majority_class)[0]
minority_indices = np.where(y == minority_class)[0]

# Downsample majority class to match the size of the minority class
majority_downsampled = resample(majority_indices, 
                                n_samples=len(minority_indices),
                                random_state=42)

# Combine indices to create a balanced dataset
balanced_indices = np.concatenate([minority_indices, majority_downsampled])

X_balanced = X.iloc[balanced_indices]
y_balanced = y[balanced_indices]

# Define a custom scoring function for specificity
def specificity_score(y_true, y_pred):
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    return tn / (tn + fp)

# Define a custom scoring function for classification error
def classification_error(y_true, y_pred):
    return np.mean(y_true != y_pred)

# Define a dictionary of scoring metrics
scoring = {
    'auc': 'roc_auc',
    'ce': make_scorer(classification_error, greater_is_better=False),
    'acc': 'accuracy',
    'precision': 'precision',
    'recall': 'recall',
    'sensitivity': 'recall',
    'specificity': make_scorer(specificity_score)
}

# Define a list of models with their respective pipelines
models = [
    ('XGBoost', Pipeline([
        ('feature_selection', SelectKBest(f_classif, k=50)),
        ('xgboost', xgb.XGBClassifier(use_label_encoder=False, eval_metric='logloss',
                                      subsample=0.8, n_estimators=1000, min_child_weight=1,
                                      max_depth=5, learning_rate=0.1, gamma=0, colsample_bytree=1.0))
    ])),
    ('RandomForest', Pipeline([
        ('feature_selection', SelectKBest(f_classif, k=10)),
        ('randomforest', RandomForestClassifier(n_estimators=100, min_samples_split=5,
                                                min_samples_leaf=1, max_features='log2',
                                                max_depth=None, class_weight='balanced'))
    ])),
    ('NaiveBayes', Pipeline([
        ('scaler', RobustScaler()),
        ('feature_selection', SelectKBest(f_classif)),
        ('pca', PCA()),
        ('nb', GaussianNB())
    ])),
    ('Logistic', Pipeline([
        ('scaler', RobustScaler()),
        ('feature_selection', RFE(estimator=LogisticRegression(), n_features_to_select=27)),
        ('logistic', LogisticRegression(C=94.88955372533333, class_weight='balanced', 
                                        l1_ratio=0.3854165025399161, penalty='l2', 
                                        solver='saga', max_iter=10000))
    ])),
    ('LDA', Pipeline([
        ('scaler', RobustScaler()),
        ('feature_selection', SelectKBest(f_classif, k=84)),
        ('pca', PCA(n_components=5)),
        ('lda', LinearDiscriminantAnalysis(solver='lsqr', shrinkage=0.3051428023145554))
    ])),
    ('KNN', Pipeline([
        ('feature_selection', SelectKBest(f_classif, k=10)),
        ('knn', KNeighborsClassifier(weights='distance', p=2, n_neighbors=5))
    ])),
    ('GradientBoosting', Pipeline([
        ('feature_selection', SelectKBest(f_classif, k=10)),
        ('gradientboosting', GradientBoostingClassifier(n_estimators=500, min_samples_split=5,
                                                        min_samples_leaf=1, max_depth=7, learning_rate=0.3))
    ])),
    ('DecisionTree', Pipeline([
        ('feature_selection', SelectKBest(f_classif, k=10)),
        ('decisiontree', DecisionTreeClassifier(min_samples_split=5, min_samples_leaf=2,
                                                max_features='log2', max_depth=10,
                                                criterion='entropy', class_weight=None))
    ])),
    ('NNET', Pipeline([
        ('scaler', StandardScaler()),
        ('feature_selection', SelectKBest(f_classif, k=50)),
        ('nnet', MLPClassifier(activation='relu', alpha=0.001,
                               hidden_layer_sizes=(50, 25),
                               learning_rate_init=0.0001,
                               max_iter=10000))
    ])),
    ('SVM', Pipeline([
        ('scaler', StandardScaler()),
        ('feature_selection', SelectKBest(f_classif, k=30)),
        ('svm', SVC(C=0.1, gamma=0.1, kernel='poly', probability=True))
    ]))
]

# Set up the outer cross-validation strategy
cv_outer = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

# Initialize lists to store results and model names
results = []
model_names = []

# Evaluate each model
for name, model in models:
    print(f"Evaluating model: {name}")
    
    # Perform cross-validation
    scores = cross_validate(model, X_balanced, y_balanced, cv=cv_outer, scoring=scoring)
    
    # Calculate mean scores for each metric
    mean_scores = {metric: np.mean(score) for metric, score in scores.items() if metric.startswith('test_')}
    results.append(list(mean_scores.values()))
    model_names.append(name)
    
    print(f"Evaluation completed for {name}")
    print("-" * 50)

# Create a DataFrame from the results
if results:
    results_df = pd.DataFrame(results, columns=[metric[5:] for metric in mean_scores.keys()], index=model_names)
    results_df['ce'] = results_df['ce'].abs()  # Take the absolute value of classification error

    # Plot a heatmap of the model performance
    plt.figure(figsize=(15, 10))
    sns.heatmap(results_df, annot=True, cmap="YlGnBu", fmt='.4f', vmin=0, vmax=1.0)
    plt.title("Model Performance Heatmap")
    plt.tight_layout()
    plt.savefig('model_performance_heatmap_inter.pdf')
    plt.close()

# Print the performance metrics for each model
print("\nPerformance for each model:")
for name, scores in zip(model_names, results):
    print(f"{name}:")
    for metric, score in zip([metric[5:] for metric in mean_scores.keys()], scores):
        if metric == 'ce':
            print(f"  {metric}: {abs(score):.4f}")
        else:
            print(f"  {metric}: {score:.4f}")
    print()