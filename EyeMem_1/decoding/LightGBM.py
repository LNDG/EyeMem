#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 16:09:42 2025

@author: kloosterman
"""

import lightgbm as lgb
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score

params = {
    'objective': 'multiclass',  # Multi-class classification
    'num_class': 150,           # Number of classes
    'metric': 'multi_logloss',  # Loss function for probability predictions
    'boosting_type': 'gbdt',    # Gradient Boosting Decision Trees
    'num_leaves': 64,           # Number of leaves in a tree (higher = complex model)
    'learning_rate': 0.05,      # Controls how fast the model learns
    'feature_fraction': 0.8,    # Uses 80% of features per tree (for regularization)
    'bagging_fraction': 0.8,    # Uses 80% of data per tree (for regularization)
    'bagging_freq': 5,          # Perform bagging every 5 iterations
    'verbose': -1               # Suppress warnings
}

# Load fMRI data (replace with actual file)
df = pd.read_csv("fmri_features.csv")

# Extract features and labels
X = df.iloc[:, 2:-1].values  # Features
y = df.iloc[:, -1].values    # Labels (150 classes)

# Train-test split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, stratify=y)

# Create LightGBM dataset
train_data = lgb.Dataset(X_train, label=y_train)
test_data = lgb.Dataset(X_test, label=y_test, reference=train_data)

# Train LightGBM model
model = lgb.train(params, train_data, valid_sets=[test_data], early_stopping_rounds=50, verbose_eval=False)

# Predict probabilities
y_pred_probs = model.predict(X_test)

# Convert probabilities to class predictions
y_pred = np.argmax(y_pred_probs, axis=1)

# Evaluate accuracy
acc = accuracy_score(y_test, y_pred)
print(f"Accuracy: {acc:.4f}")
